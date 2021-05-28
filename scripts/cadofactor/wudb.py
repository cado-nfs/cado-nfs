#!/usr/bin/env python3

# TODO:
# FILES table: OBSOLETE column
#     OBSOLETE says that this file was replaced by a newer version, for example checkrels may want to create a new file with only part of the output.
#     Should this be a pointer to the file that replaced the obsolete one? How to signal a file that is obsolete, but not replaced by anything?
#     If one file is replaced by several (say, due to a data corruption in the middle), we need a 1:n relationship. If several files are replaced by one (merge),
#     we need n:1. What do? Do we really want an n:n relationship here? Disallow fragmenting files, or maybe simply not track it in the DB if we do?
# FILES table: CHECKSUM column
#     We need a fast check that the information stored in the DB still accurately reflects the file system contents. The test should also warn about files in upload/ which are not listed in DB


# Make Python 2.7 use the print() syntax from Python 3
from __future__ import print_function

import sys
import sqlite3
try:
    import mysql
    import mysql.connector
    HAVE_MYSQL=True
except ImportError:
    HAVE_MYSQL=False

import threading
import traceback
import collections
import abc
from datetime import datetime
import re
import time

if re.search("^/", "@CMAKE_INSTALL_PREFIX@/@LIBSUFFIX@"):
    sys.path.append("@CMAKE_INSTALL_PREFIX@/@LIBSUFFIX@")

from workunit import Workunit
if sys.version_info.major == 3:
    from queue import Queue
else:
    from Queue import Queue
import patterns
import cadologger
import logging
from collections.abc import MutableMapping
from collections.abc import Mapping

DEBUG = 1
exclusive_transaction = [None, None]

DEFERRED = object()
IMMEDIATE = object()
EXCLUSIVE = object()

logger = logging.getLogger("Database")
logger.setLevel(logging.NOTSET)


PRINTED_CANCELLED_WARNING = False

def join3(l, pre=None, post=None, sep=", "):
    """
    If any parameter is None, it is interpreted as the empty string
    >>> join3 ( ('a'), pre="+", post="-", sep=", ")
    '+a-'
    >>> join3 ( ('a', 'b'), pre="+", post="-", sep=", ")
    '+a-, +b-'
    >>> join3 ( ('a', 'b'))
    'a, b'
    >>> join3 ( ('a', 'b', 'c'), pre="+", post="-", sep=", ")
    '+a-, +b-, +c-'
    """
    if pre is None:
        pre = ""
    if post is None:
        post = ""
    if sep is None:
        sep = "";
    return sep.join([pre + k + post for k in l])

def dict_join3(d, sep=None, op=None, pre=None, post=None):
    """
    If any parameter is None, it is interpreted as the empty string
    >>> t = dict_join3 ( {"a": "1", "b": "2"}, sep=",", op="=", pre="-", post="+")
    >>> t == '-a=1+,-b=2+' or t == '-b=2+,-a=1+'
    True
    """
    if pre is None:
        pre = ""
    if post is None:
        post = ""
    if sep is None:
        sep = "";
    if op is None:
        op = ""
    return sep.join([pre + op.join(k) + post for k in d.items()])

def conn_commit(conn):
    logger.transaction("Commit on connection %d", id(conn))
    if DEBUG > 1:
        if not exclusive_transaction[0] is None and not conn is exclusive_transaction[0]:
            logger.warning("Commit on connection %d, but exclusive lock was on %d", id(conn), id(exclusive_transaction[0]))
        exclusive_transaction[0] = None
        exclusive_transaction[1] = None
    conn.commit()

def conn_close(conn):
    # I'm really having difficulties here. I can't see what's going on.
    # Sometimes I have an uncaught exception popping up.
    #target = 92800609832959449330691138186
    #log(target) = 32359472153599817010011705
    #Warning:Database: Connection 140584385754280 being closed while in transaction
    #Exception ignored in: <bound method WuAccess.__del__ of <wudb.WuAccess object at 0x7fdc5a5fb470>>
    #Traceback (most recent call last):
    #  File "/home/thome/NFS/cado/scripts/cadofactor/wudb.py", line 1128, in __del__
    #  File "/home/thome/NFS/cado/scripts/cadofactor/wudb.py", line 107, in conn_close
    #  File "/usr/lib/python3.5/logging/__init__.py", line 1292, in warning
    #  File "/usr/lib/python3.5/logging/__init__.py", line 1416, in _log
    #  File "/usr/lib/python3.5/logging/__init__.py", line 1426, in handle
    #  File "/usr/lib/python3.5/logging/__init__.py", line 1488, in callHandlers
    #  File "/usr/lib/python3.5/logging/__init__.py", line 856, in handle
    #  File "/usr/lib/python3.5/logging/__init__.py", line 1048, in emit
    #  File "/usr/lib/python3.5/logging/__init__.py", line 1038, in _open
    #NameError: name 'open' is not defined
    #
    try:
        logger.transaction("Closing connection %d", id(conn))
        if conn.in_transaction:
            logger.warning("Connection %d being closed while in transaction", id(conn))
        conn.close()
    except:
        pass


# Dummy class for defining "constants" with reverse lookup
STATUS_NAMES = ["AVAILABLE", "ASSIGNED", "NEED_RESUBMIT", "RECEIVED_OK",
         "RECEIVED_ERROR", "VERIFIED_OK", "VERIFIED_ERROR", "CANCELLED"]
STATUS_VALUES = range(len(STATUS_NAMES))
WuStatusBase = collections.namedtuple("WuStatusBase", STATUS_NAMES)
class WuStatusClass(WuStatusBase):
    def check(self, status):
        assert status in self
    def get_name(self, status):
        self.check(status)
        return STATUS_NAMES[status]

WuStatus = WuStatusClass(*STATUS_VALUES)


def check_tablename(name):
    """ Test whether name is a valid SQL table name.

    Raise an exception if it isn't.
    """
    no_ = name.replace("_", "")
    if not no_[0].isalpha() or not no_[1:].isalnum():
        raise Exception("%s is not valid for an SQL table name" % name)

# If we try to update the status in any way other than progressive
# (AVAILABLE -> ASSIGNED -> ...), we raise this exception
class StatusUpdateError(Exception):
    pass

# I wish I knew how to make that inherit from a template argument (which
# would be sqlite3.Cursor or mysql.Cursor). I'm having difficulties to
# grok that syntax though, so let's stay simple and stupid. We'll have a
# *member* which is the cursor object, and so be it.
class CursorWrapperBase(object,metaclass=abc.ABCMeta):
    """ This class represents a DB cursor and provides convenience functions
        around SQL queries. In particular it is meant to provide an
        (1) an interface to SQL functionality via method calls with parameters,
        and
        (2) hiding some particularities of the SQL variant of the underlying
            DBMS as far as possible """

    # This is used in where queries; it converts from named arguments such as
    # "eq" to a binary operator such as "="
    name_to_operator = {"lt": "<", "le": "<=", "eq": "=", "ge": ">=", "gt" : ">", "ne": "!=", "like": "like"}
    @abc.abstractproperty
    def cursor(self):
        pass

    @abc.abstractproperty
    def connection(self):
        pass

    # override in the derived cursor class if needed
    @property
    def _string_translations(self):
        return []

    # override in the derived cursor class if needed
    def translations(self, x):
        if type(x) == tuple:
            return tuple([self.translations(u) for u in x])
        elif type(x) == list:
            return [self.translations(u) for u in x]
        else:
            v=x
            for a,b in self._string_translations:
                v,nrepl=re.subn(a, b, v)
            return v

    # override in the derived cursor class if needed
    @property
    def parameter_auto_increment(self):
        return "?"

    def __init__(self):
        pass

    def in_transaction(self):
        return self.connection.in_transaction

    @staticmethod
    def _without_None(d):
        """ Return a copy of the dictionary d, but without entries whose values
            are None """
        return {k[0]:k[1] for k in d.items() if k[1] is not None}

    @staticmethod
    def as_string(d):
        if d is None:
            return ""
        else:
            return ", " + dict_join3(d, sep=", ", op=" AS ")

    def _where_str(self, name, **args):
        where = ""
        values = []
        qm=self.parameter_auto_increment
        for opname in args:
            if args[opname] is None:
                continue
            if where == "":
                where = " " + name + " "
            else:
                where = where + " AND "
            where = where + join3(args[opname].keys(),
                        post=" " + self.name_to_operator[opname] + " " + qm,
                        sep=" AND ")
            values = values + list(args[opname].values())
        return (where, values)

    def _exec(self, command, values=None):
        """ Wrapper around self.execute() that prints arguments
            for debugging and retries in case of "database locked" exception """

        # FIXME: should be the caller's class name, as _exec could be
        # called from outside this class
        classname = self.__class__.__name__
        parent = sys._getframe(1).f_code.co_name
        command = self.translations(command)
        command_str = command.replace("?", "%r")
        if not values is None:
            command_str = command_str % tuple(values)
        logger.transaction("%s.%s(): connection = %s, command = %s",
                           classname, parent, id(self.connection), command_str)
        i = 0
        while True:
            try:
                if values is None or len(values)==0:
                    self.cursor.execute(command)
                else:
                    self.cursor.execute(command, values)
                break
            except (sqlite3.OperationalError, sqlite3.DatabaseError) as e:
                if str(e) == "database disk image is malformed" or \
                        str(e) == "disk I/O error":
                    logger.critical("sqlite3 reports error accessing the database.")
                    logger.critical("Database file may have gotten corrupted, "
                            "or maybe filesystem does not properly support "
                            "file locking.")
                    raise
                if str(e) != "database is locked":
                   raise
                i += 1
                time.sleep(1) # wait for 1 second if database is locked
                if i == 10:
                    logger.critical("You might try 'fuser xxx.db' to see which process is locking the database")
                    raise
        logger.transaction("%s.%s(): connection = %s, command finished",
                           classname, parent, id(self.connection))

    def begin(self, mode=None):
        if mode is None:
            self._exec("BEGIN")
        elif mode is DEFERRED:
            self._exec("BEGIN DEFERRED")
        elif mode is IMMEDIATE:
            self._exec("BEGIN IMMEDIATE")
        elif mode is EXCLUSIVE:
            if DEBUG > 1:
                tb = traceback.extract_stack()
                if not exclusive_transaction == [None, None]:
                    old_tb_str = "".join(traceback.format_list(exclusive_transaction[1]))
                    new_tb_str = "".join(traceback.format_list(tb))
                    logger.warning("Called Cursor.begin(EXCLUSIVE) when there was aleady an exclusive transaction %d\n%s",
                                        id(exclusive_transaction[0]), old_tb_str)
                    logger.warning("New transaction: %d\n%s", id(self.connection), new_tb_str)

            self._exec("BEGIN EXCLUSIVE")

            if DEBUG > 1:
                assert exclusive_transaction == [None, None]
                exclusive_transaction[0] = self.connection
                exclusive_transaction[1] = tb
        else:
            raise TypeError("Invalid mode parameter: %r" % mode)

    def pragma(self, prag):
        self._exec("PRAGMA %s;" % prag)

    def create_table(self, table, layout):
        """ Creates a table with fields as described in the layout parameter """
        command = "CREATE TABLE IF NOT EXISTS %s( %s );" % \
                  (table, ", ".join(" ".join(k) for k in layout))
        self._exec (command)

    def create_index(self, name, table, columns):
        # we get so many of these...
        try:
            """ Creates an index with fields as described in the columns list """
            command = self.translations("CREATE INDEX IF NOT EXISTS") + " %s ON %s( %s );" % (name, table, ", ".join(columns))
            self._exec (command)
        except Exception as e:
            logger.warning(e)
            pass

    def insert(self, table, d):
        """ Insert a new entry, where d is a dictionary containing the
            field:value pairs. Returns the row id of the newly created entry """
        # INSERT INTO table (field_1, field_2, ..., field_n)
        # 	VALUES (value_1, value_2, ..., value_n)

        # Fields is a copy of d but with entries removed that have value None.
        # This is done primarily to avoid having "id" listed explicitly in the
        # INSERT statement, because the DB fills in a new value automatically
        # if "id" is the primary key. But I guess not listing field:NULL items
        # explicitly in an INSERT is a good thing in general
        fields = self._without_None(d)
        fields_str = ", ".join(fields.keys())

        qm=self.parameter_auto_increment
        sqlformat = ", ".join((qm,) * len(fields)) # sqlformat = "?, ?, ?, " ... "?"
        command = "INSERT INTO %s( %s ) VALUES ( %s );" \
                  % (table, fields_str, sqlformat)
        values = list(fields.values())
        self._exec(command, values)
        rowid = self.lastrowid
        return rowid

    def update(self, table, d, **conditions):
        """ Update fields of an existing entry. conditions specifies the where
            clause to use for to update, entries in the dictionary d are the
            fields and their values to update """
        # UPDATE table SET column_1=value1, column2=value_2, ...,
        # column_n=value_n WHERE column_n+1=value_n+1, ...,
        qm=self.parameter_auto_increment
        setstr = join3(d.keys(), post = " = " + qm, sep = ", ")
        (wherestr, wherevalues) = self._where_str("WHERE", **conditions)
        command = "UPDATE %s SET %s %s" % (table, setstr, wherestr)
        values = list(d.values()) + wherevalues
        self._exec(command, values)

    def where_query(self, joinsource, col_alias=None, limit=None, order=None,
                    **conditions):
        # Table/Column names cannot be substituted, so include in query directly.
        (WHERE, values) = self._where_str("WHERE", **conditions)
        if order is None:
            ORDER = ""
        else:
            if not order[1] in ("ASC", "DESC"):
                raise Exception
            ORDER = " ORDER BY %s %s" % (order[0], order[1])
        if limit is None:
            LIMIT = ""
        else:
            LIMIT = " LIMIT %s" % int(limit)
        AS = self.as_string(col_alias);
        command = "SELECT * %s FROM %s %s %s %s" \
                  % (AS, joinsource, WHERE, ORDER, LIMIT)
        return (command, values)

    def where(self, joinsource, col_alias=None, limit=None, order=None,
              values=[], **conditions):
        """ Get a up to "limit" table rows (limit == 0: no limit) where
            the key:value pairs of the dictionary "conditions" are set to the
            same value in the database table """
        (command, newvalues) = self.where_query(joinsource, col_alias, limit,
                                             order, **conditions)
        self._exec(command + ";", values + newvalues)

    def count(self, joinsource, **conditions):
        """ Count rows where the key:value pairs of the dictionary "conditions" are
            set to the same value in the database table """

        # Table/Column names cannot be substituted, so include in query directly.
        (WHERE, values) = self._where_str("WHERE", **conditions)

        command = "SELECT COUNT(*) FROM %s %s;" % (joinsource, WHERE)
        self._exec(command, values)
        r = self.cursor.fetchone()
        return int(r[0])

    def delete(self, table, **conditions):
        """ Delete the rows specified by conditions """
        (WHERE, values) = self._where_str("WHERE", **conditions)
        command = "DELETE FROM %s %s;" % (table, WHERE)
        self._exec(command, values)

    def where_as_dict(self, joinsource, col_alias=None, limit=None,
                      order=None, values=[], **conditions):
        self.where(joinsource, col_alias=col_alias, limit=limit,
                      order=order, values=values, **conditions)
        # cursor.description is a list of lists, where the first element of
        # each inner list is the column name
        result = []
        desc = [k[0] for k in self.cursor.description]
        row = self.cursor.fetchone()
        while row is not None:
            # print("Cursor.where_as_dict(): row = %s" % row)
            result.append(dict(zip(desc, row)))
            row = self.cursor.fetchone()
        return result
    def execute(self, *args, **kwargs):
        return self._exec(*args, **kwargs)
    def fetchone(self, *args, **kwargs):
        return self.cursor.fetchone(*args, **kwargs)
    def close(self):
        self.cursor.close()
    @property
    def lastrowid(self):
        self.cursor.lastrowid

class DB_base(object):
    @property
    def general_pattern(self):
        return r"(?:db:)?(\w+)://(?:(?:(\w+)(?::(.*))?@)?(?:([\w\.]+)|\[([\d:]+)*\])(?::(\d+))?/)?(.*)$"
    def __init__(self, uri, backend_pattern=None):
        self.uri = uri
        foo=re.match(self.general_pattern,uri)
        if not foo:
            raise ValueError("db URI %s does not match regexp %s" % (uri,self.general_pattern))
        self.hostname=foo.group(4)
        self.host_ipv6=False
        if not self.hostname:
            self.hostname=foo.group(5)
            self.host_ipv6=True
        self.backend=foo.group(1)
        if backend_pattern is not None and not re.match(backend_pattern, self.backend):
            raise ValueError("back-end type %s not supported, expected %s" % (self.backend, backend_pattern))
        self.db_connect_args=dict(
                user=foo.group(2),
                password=foo.group(3),
                host=self.hostname,
                port=foo.group(6)
        )
        self.db_name=foo.group(7)
        self.talked=False
        # logger.info("Database URI is %s" % self.uri_without_credentials)
    @property
    def uri_without_credentials(self):
        text="db:%s://" % self.backend
        d=self.db_connect_args
        if "host" in d:
            if "user" in d:
                text+="USERNAME"
                if "password" in d:
                    text+=":PASSWORD"
                text+="@"
            if self.host_ipv6:
                text+="[%s]" % d["host"]
            else:
                text+=d["host"]
            if "port" in d:
                text+=":%s" % d["port"]
            text+="/"
        text+=self.db_name
        return text
    def advertise_connection(self):
        if not self.talked:
            logger.info("Opened connection to database %s" % self.db_name)
            self.talked=True

class DB_SQLite(DB_base):
    class CursorWrapper(CursorWrapperBase):
        @property
        def cursor(self):
            return self.__cursor
        @property
        def connection(self):
            return self.cursor.connection
        def __init__(self, cursor, *args, **kwargs):
            self.__cursor=cursor
            super().__init__(*args, **kwargs)
    class ConnectionWrapper(sqlite3.Connection):
        def cursor(self):
            return DB_SQLite.CursorWrapper(super().cursor())
        def __init__(self, *args, **kwargs):
            super().__init__(isolation_level=None, *args, **kwargs)
    def connect(self):
        c=self.ConnectionWrapper(self.path)
        self.advertise_connection()
        return c
    # FIXME I think that in any case the sqlite3 module ends up creating
    # the db, no ?
    def __init__(self, uri, create=False):
        super().__init__(uri, backend_pattern="sqlite3?")
        self.path = self.db_name

if HAVE_MYSQL:
    class DB_MySQL(DB_base):
        class CursorWrapper(CursorWrapperBase):
            @property
            def parameter_auto_increment(self):
                return "%s"
            @property
            def _string_translations(self):
                return [
                        ('\\bASC\\b', "AUTO_INCREMENT"),
                        ('\\bCREATE INDEX IF NOT EXISTS\\b', "CREATE INDEX"),
                        ('\\bBEGIN EXCLUSIVE\\b', "START TRANSACTION"),
                        ('\\bpurge\\b', "purgetable"),
                ]
            @property
            def cursor(self):
                return self.__cursor
            @property
            def connection(self):
                return self._connection
            def __init__(self, cursor, connection=None, *args, **kwargs):
                self._connection = connection
                self.__cursor=cursor
                super().__init__(*args, **kwargs)

        class ConnectionWrapper(object):
            def _reconnect_anonymous(self):
                self._conn = mysql.connector.connect(**self._db_factory.db_connect_args)
            def _reconnect(self):
                self._conn = mysql.connector.connect(database=self._db_factory.db_name, **self._db_factory.db_connect_args)
            def cursor(self):
                # provide some retry capability. This must be done on the
                # connection object, since reconnecting changes the
                # connection member.
                for i in range(10):
                    try:
                        c = self._conn.cursor()
                        break
                    except mysql.connector.errors.OperationalError as e:
                        logger.warning("Got exception connecting to the database, retrying (#%d)" % i)
                        if self.db:
                            self._reconnect()
                        else:
                            raise
                self._conn.commit()
                return DB_MySQL.CursorWrapper(c, connection=self)
            def __init__(self, db_factory, create=False):
                self._db_factory = db_factory
                db_name = self._db_factory.db_name
                if create:
                    try:
                        self._reconnect()
                    except mysql.connector.errors.ProgrammingError:
                        # need to create the database first. Do it by
                        # hand, with a connection which starts without a
                        # database name.
                        logger.info("Creating database %s" % db_name)
                        self._reconnect_anonymous()
                        cursor = self._conn.cursor()
                        cursor.execute("CREATE DATABASE %s;" % db_name)
                        cursor.execute("USE %s;" % db_name)
                        cursor.execute("SET autocommit = 1")
                        self._conn.commit()
                else:
                    self._reconnect()
            def rollback(self):
                self._conn.rollback()
            def close(self):
                self._conn.close()
            def commit(self):
                self._conn.commit()
            @property
            def in_transaction(self):
                return self._conn.in_transaction
        def connect(self, *args, **kwargs):
            return self.ConnectionWrapper(self, *args, **kwargs)
        def __init__(self, uri,create=False):
            super().__init__(uri, backend_pattern="mysql")
            self.path = None
            if create:
                conn=self.connect(create=True)
                conn.close()

class DBFactory(object):
    # This class initializes the database from the supplied db uri.
    # db:engine:[//[user[:password]@][host][:port]/][dbname][?params][#fragment]
    def __init__(self, uri, *args, **kwargs):
        self.uri = uri
        self.base = None
        error={}
        sc=DB_base.__subclasses__()
        for c in sc:
            # logger.info("Trying database back-end %s (among %d)" % (c, len(sc)))
            try:
                self.base = c(uri, *args, **kwargs)
                break
            except ValueError as err:
                error[str(c)]=err
                pass
        if self.base is None:
            msg = "Cannot use database URI %s" % uri
            msg += "\n" + "Messages received from %d backends:" % len(sc)
            for c in error.keys():
                msg += "\n" + "Error from %s: %s" % (c, error[c])
            raise ValueError(msg)
    def connect(self):
        return self.base.connect()
    @property
    def uri_without_credentials(self):
        return self.base.uri_without_credentials
    @property
    def path(self):
        # TODO: remove
        return self.base.path


class DbTable(object):
    """ A class template defining access methods to a database table """

    @staticmethod
    def _subdict(d, l):
        """ Returns a dictionary of those key:value pairs of d for which key
            exists l """
        if d is None:
            return None
        return {k:d[k] for k in d.keys() if k in l}

    def _get_colnames(self):
        return [k[0] for k in self.fields]

    def getname(self):
        return self.tablename

    def getpk(self):
        return self.primarykey

    def dictextract(self, d):
        """ Return a dictionary with all those key:value pairs of d
            for which key is in self._get_colnames() """
        return self._subdict(d, self._get_colnames())

    def create(self, cursor):
        fields = list(self.fields)
        if self.references:
            # If this table references another table, we use the primary
            # key of the referenced table as the foreign key name
            r = self.references # referenced table
            fk = (r.getpk(), "INTEGER", "REFERENCES %s ( %s ) " \
                  % (r.getname(), r.getpk()))
            fields.append(fk)
        cursor.create_table(self.tablename, fields)
        if self.references:
            # We always create an index on the foreign key
            cursor.create_index(self.tablename + "_pkindex", self.tablename,
                                (fk[0], ))
        for indexname in self.index:
            # cursor.create_index(self.tablename + "_" + indexname, self.tablename, self.index[indexname])
            try:
                cursor.create_index(self.tablename + "_" + indexname + "_index",
                                    self.tablename, self.index[indexname])
            except Exception as e:
                logger.warning(e)
                pass

    def insert(self, cursor, values, foreign=None):
        """ Insert a new row into this table. The column:value pairs are
            specified key:value pairs of the dictionary d.
            The database's row id for the new entry is stored in
            d[primarykey] """
        d = self.dictextract(values)
        assert self.primarykey not in d or d[self.primarykey] is None
        # If a foreign key is specified in foreign, add it to the column
        # that is marked as being a foreign key
        if foreign:
            r = self.references.primarykey
            assert not r in d or d[r] is None
            d[r] = foreign
        values[self.primarykey] = cursor.insert(self.tablename, d)

    def insert_list(self, cursor, values, foreign=None):
        for v in values:
            self.insert(cursor, v, foreign)

    def update(self, cursor, d, **conditions):
        """ Update an existing row in this table. The column:value pairs to
            be written are specified key:value pairs of the dictionary d """
        cursor.update(self.tablename, d, **conditions)

    def delete(self, cursor, **conditions):
        """ Delete an existing row in this table """
        cursor.delete(self.tablename, **conditions)

    def where(self, cursor, limit=None, order=None, **conditions):
        assert order is None or order[0] in self._get_colnames()
        return cursor.where_as_dict(self.tablename, limit=limit,
                                    order=order, **conditions)


class WuTable(DbTable):
    tablename = "workunits"
    fields = (
        ("wurowid", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"),
        ("wuid", "VARCHAR(512)", "UNIQUE NOT NULL"),
        ("submitter", "VARCHAR(512)", ""),
        ("status", "INTEGER", "NOT NULL"),
        ("wu", "TEXT", "NOT NULL"),
        ("timecreated", "TEXT", ""),
        ("timeassigned", "TEXT", ""),
        ("assignedclient", "TEXT", ""),
        ("timeresult", "TEXT", ""),
        ("resultclient", "TEXT", ""),
        ("errorcode", "INTEGER", ""),
        ("failedcommand", "INTEGER", ""),
        ("timeverified", "TEXT", ""),
        ("retryof", "INTEGER", "REFERENCES %s" % tablename),
        ("priority", "INTEGER", "")
    )
    primarykey = fields[0][0]
    references = None
    index = {"wuid": (fields[1][0],),
             "submitter" : (fields[2][0],),
             "priority" : (fields[14][0],),
             "status" : (fields[3][0],)
    }

class FilesTable(DbTable):
    tablename = "files"
    fields = (
        ("filesrowid", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"),
        ("filename", "TEXT", ""),
        ("path", "VARCHAR(512)", "UNIQUE NOT NULL"),
        ("type", "TEXT", ""),
        ("command", "INTEGER", "")
    )
    primarykey = fields[0][0]
    references = WuTable()
    index = {}


# The sqrt_factors table contains the input number to be factored. As
# such, we must make sure that it's permitted to go at least as far as we
# intend to go. 200 digits is definitely too small.
class DictDbTable(DbTable):
    fields = (
        ("rowid", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"),
        ("kkey", "VARCHAR(300)", "UNIQUE NOT NULL"),
        ("type", "INTEGER", "NOT NULL"),
        ("value", "TEXT", "")
        )
    primarykey = fields[0][0]
    references = None
    def __init__(self, *args, name = None, **kwargs):
        self.tablename = name
        # index creation now always prepends the table name, and appends "index"
        self.index = {"dictdb_kkey": ("kkey",)} # useful ?
        super().__init__(*args, **kwargs)


class DictDbAccess(MutableMapping):
    """ A DB-backed flat dictionary.

    Flat means that the value of each dictionary entry must be a type that
    the underlying DB understands, like integers, strings, etc., but not
    collections or other complex types.

    A copy of all the data in the table is kept in memory; read accesses
    are always served from the in-memory dict. Write accesses write through
    to the DB.

    >>> conn = DBFactory('db:sqlite3://:memory:').connect()
    >>> d = DictDbAccess(conn, 'test')
    >>> d == {}
    True
    >>> d['a'] = '1'
    >>> d == {'a': '1'}
    True
    >>> d['a'] = 2
    >>> d == {'a': 2}
    True
    >>> d['b'] = '3'
    >>> d == {'a': 2, 'b': '3'}
    True
    >>> del(d)
    >>> d = DictDbAccess(conn, 'test')
    >>> d == {'a': 2, 'b': '3'}
    True
    >>> del(d['b'])
    >>> d == {'a': 2}
    True
    >>> d.setdefault('a', '3')
    2
    >>> d == {'a': 2}
    True
    >>> d.setdefault('b', 3.0)
    3.0
    >>> d == {'a': 2, 'b': 3.0}
    True
    >>> d.setdefault(None, {'a': '3', 'c': '4'})
    >>> d == {'a': 2, 'b': 3.0, 'c': '4'}
    True
    >>> d.update({'a': '3', 'd': True})
    >>> d == {'a': '3', 'b': 3.0, 'c': '4', 'd': True}
    True
    >>> del(d)
    >>> d = DictDbAccess(conn, 'test')
    >>> d == {'a': '3', 'b': 3.0, 'c': '4', 'd': True}
    True
    >>> d.clear(['a', 'd'])
    >>> d == {'b': 3.0, 'c': '4'}
    True
    >>> del(d)
    >>> d = DictDbAccess(conn, 'test')
    >>> d == {'b': 3.0, 'c': '4'}
    True
    >>> d.clear()
    >>> d == {}
    True
    >>> del(d)
    >>> d = DictDbAccess(conn, 'test')
    >>> d == {}
    True
    """

    types = (str, int, float, bool)

    def __init__(self, db, name):
        ''' Attaches to a DB table and reads values stored therein.

        db can be a string giving the file name for the DB (same as for
        sqlite3.connect()), or an open DB connection. The latter is allowed
        primarily for making the doctest work, so we can reuse the same
        memory-backed DB connection, but it may be useful in other contexts.
        '''

        if isinstance(db, DBFactory):
            self._db = db
            self._conn = db.connect()
            self._ownconn = True
        elif isinstance(db, str):
            raise ValueError("unexpected: %s" % db)
        else:
            self._db = None
            self._conn = db
            self._ownconn = False
        self._table = DictDbTable(name = name)
        # Create an empty table if none exists
        cursor = self.get_cursor()
        self._table.create(cursor);
        # Get the entries currently stored in the DB
        self._data = self._getall()
        cursor.close()

    def get_cursor(self):
        return self._conn.cursor()

    # Implement the abstract methods defined by collections.abc.MutableMapping
    # All but __del__ and __setitem__ are simply passed through to the self._data
    # dictionary

    def __getitem__(self, key):
        return self._data.__getitem__(key)

    def __iter__(self):
        return self._data.__iter__()

    def  __len__(self):
        return self._data.__len__()

    def __str__(self):
        return self._data.__str__()

    def __del__(self):
        """ Close the DB connection and delete the in-memory dictionary """
        if self._ownconn:
            # When we shut down Python hard, e.g., in an exception, the
            # conn_close() function object may have been destroyed already
            # and trying to call it would raise another exception.
            if callable(conn_close):
                conn_close(self._conn)
            else:
                self._conn.close()

    def __convert_value(self, row):
        valuestr = row["value"]
        valuetype = row["type"]
        # Look up constructor for this type
        typecon = self.types[int(valuetype)]
        # Bool is handled separately as bool("False") == True
        if typecon == bool:
            if valuestr == "True":
                return True
            elif valuestr == "False":
                return False
            else:
                raise ValueError("Value %s invalid for Bool type", valuestr)
        return typecon(valuestr)

    def __get_type_idx(self, value):
        valuetype = type(value)
        for (idx, t) in enumerate(self.types):
            if valuetype == t:
                return idx
        raise TypeError("Type %s not supported" % str(valuetype))

    def _getall(self):
        """ Reads the whole table and returns it as a dict """
        cursor = self.get_cursor()
        rows = self._table.where(cursor)
        cursor.close()
        return {r["kkey"]: self.__convert_value(r) for r in rows}

    def __setitem_nocommit(self, cursor, key, value):
        """ Set dictionary key to value and update/insert into table,
        but don't commit. Cursor must be given
        """
        update = {"value": str(value), "type": self.__get_type_idx(value)}
        if key in self._data:
            # Update the table row where column "key" equals key
            self._table.update(cursor, update, eq={"kkey": key})
        else:
            # Insert a new row
            update["kkey"] = key
            self._table.insert(cursor, update)
        # Update the in-memory dict
        self._data[key] = value

    def __setitem__(self, key, value):
        """ Access by indexing, e.g., d["foo"]. Always commits """
        cursor = self.get_cursor()

        if not cursor.in_transaction:
            cursor.begin(EXCLUSIVE)
        self.__setitem_nocommit(cursor, key, value)
        conn_commit(self._conn)

        cursor.close()


    def __delitem__(self, key, commit=True):
        """ Delete a key from the dictionary """
        cursor = self.get_cursor()

        if not cursor.in_transaction:
            cursor.begin(EXCLUSIVE)
        self._table.delete(cursor, eq={"kkey": key})
        if commit:
            conn_commit(self._conn)

        cursor.close()
        del(self._data[key])

    def setdefault(self, key, default = None, commit=True):
        ''' Setdefault function that allows a mapping as input

        Values from default dict are merged into self, *not* overwriting
        existing values in self '''
        if key is None and isinstance(default, Mapping):
            update = {key:default[key] for key in default if not key in self}
            if update:
                self.update(update, commit=commit)
            return None
        elif not key in self:
            self.update({key:default}, commit=commit)
        return self._data[key]

    def update(self, other, commit=True):
        cursor = self.get_cursor()
        if not self._conn.in_transaction:
            cursor.begin(EXCLUSIVE)
        for (key, value) in other.items():
            self.__setitem_nocommit(cursor, key, value)
        if commit:
            conn_commit(self._conn)

        cursor.close()

    def clear(self, args = None, commit=True):
        """ Overridden clear that allows removing several keys atomically """
        cursor = self.get_cursor()
        if not self._conn.in_transaction:
            cursor.begin(EXCLUSIVE)
        if args is None:
            self._data.clear()
            self._table.delete(cursor)
        else:
            for key in args:
                del(self._data[key])
                self._table.delete(cursor, eq={"kkey": key})
        if commit:
            conn_commit(self._conn)

        cursor.close()


class Mapper(object):
    """ This class translates between application objects, i.e., Python
        directories, and the relational data layout in an SQL DB, i.e.,
        one or more tables which possibly have foreign key relationships
        that map to hierarchical data structures. For now, only one
        foreign key / subdirectory."""

    def __init__(self, table, subtables = None):
        self.table = table
        self.subtables = {}
        if subtables:
            for s in subtables.keys():
                self.subtables[s] = Mapper(subtables[s])

    def __sub_dict(self, d):
        """ For each key "k" that has a subtable assigned in "self.subtables",
        pop the entry with key "k" from "d", and store it in a new directory
        which is returned. I.e., the directory d is separated into
        two parts: the part which corresponds to subtables and is the return
        value, and the rest which is left in the input dictionary. """
        sub_dict = {}
        for s in self.subtables.keys():
            # Don't store s:None entries even if they exist in d
            t = d.pop(s, None)
            if not t is None:
                sub_dict[s] = t
        return sub_dict

    def getname(self):
        return self.table.getname()

    def getpk(self):
        return self.table.getpk()

    def create(self, cursor):
        self.table.create(cursor)
        for t in self.subtables.values():
            t.create(cursor)

    def insert(self, cursor, wus, foreign=None):
        pk = self.getpk()
        for wu in wus:
            # Make copy so sub_dict does not change caller's data
            wuc = wu.copy()
            # Split off entries that refer to subtables
            sub_dict = self.__sub_dict(wuc)
            # We add the entries in wuc only if it does not have a primary
            # key yet. If it does have a primary key, we add only the data
            # for the subtables
            if not pk in wuc:
                self.table.insert(cursor, wuc, foreign=foreign)
                # Copy primary key into caller's data
                wu[pk] = wuc[pk]
            for subtable_name in sub_dict.keys():
                self.subtables[subtable_name].insert(
                    cursor, sub_dict[subtable_name], foreign=wu[pk])

    def update(self, cursor, wus):
        pk = self.getpk()
        for wu in wus:
            assert not wu[pk] is None
            wuc = wu.copy()
            sub_dict = self.__sub_dict(wuc)
            rowid = wuc.pop(pk, None)
            if rowid:
                self.table.update(cursor, wuc, {wp: rowid})
            for s in sub.keys:
                self.subtables[s].update(cursor, sub_dict[s])

    def count(self, cursor, **cond):
        joinsource = self.table.tablename
        return cursor.count(joinsource, **cond)

    def where(self, cursor, limit = None, order = None, **cond):
        # We want:
        # SELECT * FROM (SELECT * from workunits WHERE status = 2 LIMIT 1) LEFT JOIN files USING ( wurowid );
        pk = self.getpk()
        (command, values) = cursor.where_query(self.table.tablename,
                                               limit=limit, **cond)
        joinsource = "( %s )" % command
        for s in self.subtables.keys():
            # FIXME: this probably breaks with more than 2 tables
            joinsource = "%s tmp LEFT JOIN %s USING ( %s )" \
                         % (joinsource, self.subtables[s].getname(), pk)
        # FIXME: don't get result rows as dict! Leave as tuple and
        # take them apart positionally

        rows = cursor.where_as_dict(joinsource, order=order, values=values)
        wus = []
        for r in rows:

            # Collapse rows with identical primary key
            if len(wus) == 0 or r[pk] != wus[-1][pk]:
                wus.append(self.table.dictextract(r))
                for s in self.subtables.keys():
                    wus[-1][s] = None

            for (sn, sm) in self.subtables.items():
                spk = sm.getpk()
                # if there was a match on the subtable
                if spk in r and not r[spk] is None:
                    if wus[-1][sn] == None:
                        # If this sub-array is empty, init it
                        wus[-1][sn] = [sm.table.dictextract(r)]
                    elif r[spk] != wus[-1][sn][-1][spk]:
                        # If not empty, and primary key of sub-table is not
                        # same as in previous entry, add it
                        wus[-1][sn].append(sm.table.dictextract(r))
        return wus

class WuAccess(object): # {
    """ This class maps between the WORKUNIT and FILES tables
        and a dictionary
        {"wuid": string, ..., "timeverified": string, "files": list}
        where list is None or a list of dictionaries of the from
        {"id": int, "type": int, "wuid": string, "filename": string,
        "path": string}
        Operations on instances of WuAcccess are directly carried
        out on the database persistent storage, i.e., they behave kind
        of as if the WuAccess instance were itself a persistent
        storage device """

    def __init__(self, db):
        if isinstance(db, DBFactory):
            self.conn = db.connect()
            self._ownconn = True
        elif isinstance(db, str):
            raise ValueError("unexpected")
        else:
            self.conn = db
            self._ownconn = False
        cursor = self.get_cursor()
        if isinstance(cursor, DB_SQLite.CursorWrapper):
            cursor.pragma("foreign_keys = ON")
        # I'm not sure it's relevant to do commit() at this point.
        # self.commit()
        cursor.close()
        self.mapper = Mapper(WuTable(), {"files": FilesTable()})

    def get_cursor(self):
        c = self.conn.cursor()
        return c

    def __del__(self):
        if self._ownconn:
            if callable(conn_close):
                conn_close(self.conn)
            else:
                self.conn.close()

    @staticmethod
    def to_str(wus):
        r = []
        for wu in wus:
            s = "Workunit %s:\n" % wu["wuid"]
            for (k,v) in wu.items():
                if k != "wuid" and k != "files":
                    s += "  %s: %r\n" % (k, v)
            if "files" in wu:
                s += "  Associated files:\n"
                if wu["files"] is None:
                    s += "    None\n"
                else:
                    for f in wu["files"]:
                        s += "    %s\n" % f
            r.append(s)
        return '\n'.join(r)

    @staticmethod
    def _checkstatus(wu, status):
        #logger.debug("WuAccess._checkstatus(%s, %s)", wu, status)
        wu_status = wu["status"]
        if isinstance(status, collections.abc.Container):
            ok = wu_status in status
        else:
            ok = wu_status == status
        if not ok:
            msg = "Workunit %s has status %s (%s), expected %s (%s)" % \
                  (wu["wuid"], wu_status, WuStatus.get_name(wu_status),
                   status, WuStatus.get_name(status))
            if status is WuStatus.ASSIGNED and wu_status is WuStatus.CANCELLED:
                logger.warning ("WuAccess._checkstatus(): %s, presumably timed out", msg)
                raise StatusUpdateError(msg)
            elif status is WuStatus.ASSIGNED and wu_status is WuStatus.NEED_RESUBMIT:
                logger.warning ("WuAccess._checkstatus(): %s, manually expired", msg)
                raise StatusUpdateError(msg)
            else:
                logger.error ("WuAccess._checkstatus(): %s", msg)
                raise StatusUpdateError(msg)

    # Which fields should be None for which status
    should_be_unset = {
        "errorcode": (WuStatus.AVAILABLE, WuStatus.ASSIGNED),
        "timeresult": (WuStatus.AVAILABLE, WuStatus.ASSIGNED),
        "resultclient": (WuStatus.AVAILABLE, WuStatus.ASSIGNED),
        "timeassigned": (WuStatus.AVAILABLE,),
        "assignedclient": (WuStatus.AVAILABLE,),
    }
    def check(self, data):
        status = data["status"]
        WuStatus.check(status)
        wu = Workunit(data["wu"])
        assert wu.get_id() == data["wuid"]
        if status == WuStatus.RECEIVED_ERROR:
            assert data["errorcode"] != 0
        if status == WuStatus.RECEIVED_OK:
            assert data["errorcode"] is None or data["errorcode"] == 0
        for field in self.should_be_unset:
            if status in self.should_be_unset[field]:
                assert data[field] is None

    # Here come the application-visible functions that implement the
    # "business logic": creating a new workunit from the text of a WU file,
    # assigning it to a client, receiving a result for the WU, marking it as
    # verified, or marking it as cancelled

    def _add_files(self, cursor, files, wuid=None, rowid=None):
        # Exactly one must be given
        assert not wuid is None or not rowid is None
        assert wuid is None or rowid is None
        # FIXME: allow selecting row to update directly via wuid, without
        # doing query for rowid first
        pk = self.mapper.getpk()
        if rowid is None:
            wu = get_by_wuid(cursor, wuid)
            if wu:
                rowid = wu[pk]
            else:
                return False
        colnames = ("filename", "path", "type", "command")
        # zipped length is that of shortest list, so "command" is optional
        d = (dict(zip(colnames, f)) for f in files)
        # These two should behave identically
        if True:
            self.mapper.insert(cursor, [{pk:rowid, "files": d},])
        else:
            self.mapper.subtables["files"].insert(cursor, d, foreign=rowid)

    def commit(self, do_commit=True):
        if do_commit:
            conn_commit(self.conn)

    def create_tables(self):
        cursor = self.get_cursor()
        if isinstance(cursor, DB_SQLite.CursorWrapper):
            cursor.pragma("journal_mode=WAL")
        self.mapper.create(cursor)
        self.commit()
        cursor.close()

    def _create1(self, cursor, wutext, priority=None):
        d = {
            "wuid": Workunit(wutext).get_id(),
            "wu": wutext,
            "status": WuStatus.AVAILABLE,
            "timecreated": str(datetime.utcnow())
            }
        if not priority is None:
            d["priority"] = priority
        # Insert directly into wu table
        self.mapper.table.insert(cursor, d)

    def create(self, wus, priority=None, commit=True):
        """ Create new workunits from wus which contains the texts of the
            workunit files """
        cursor = self.get_cursor()
        # todo restore transactions
        if not self.conn.in_transaction:
            cursor.begin(EXCLUSIVE)
        if isinstance(wus, str):
            self._create1(cursor, wus, priority)
        else:
            for wu in wus:
                self._create1(cursor, wu, priority)
        self.commit(commit)
        cursor.close()

    def assign(self, clientid, commit=True, timeout_hint=None):
        """ Finds an available workunit and assigns it to clientid.
            Returns the text of the workunit, or None if no available
            workunit exists """
        cursor = self.get_cursor()
        if not self.conn.in_transaction:
            cursor.begin(EXCLUSIVE)
# This "priority" stuff is the root cause for the server taking time to
# hand out WUs when the count of available WUs drops to zero.
# (introduced in 90ae4beb7 -- it's an optional-and-never-used feature
# anyway)
#        r = self.mapper.table.where(cursor, limit = 1,
#                                    order=("priority", "DESC"),
#                                    eq={"status": WuStatus.AVAILABLE})
        r = self.mapper.table.where(cursor, limit = 1,
                                    eq={"status": WuStatus.AVAILABLE})
        assert len(r) <= 1
        if len(r) == 1:
            try:
                self._checkstatus(r[0], WuStatus.AVAILABLE)
            except StatusUpdateError:
                self.commit(commit)
                cursor.close()
                raise
            if DEBUG > 0:
                self.check(r[0])
            d = {"status": WuStatus.ASSIGNED,
                 "assignedclient": clientid,
                 "timeassigned": str(datetime.utcnow())
                 }
            pk = self.mapper.getpk()
            self.mapper.table.update(cursor, d, eq={pk:r[0][pk]})
            result = r[0]["wu"]
            if timeout_hint:
                dltext = "%d\n" % int(time.time() + int(timeout_hint))
                result = result + "DEADLINE " + dltext

        else:
            result = None

        self.commit(commit)

        cursor.close()
        return result

    def get_by_wuid(self, cursor, wuid):
        r = self.mapper.where(cursor, eq={"wuid": wuid})
        assert len(r) <= 1
        if len(r) == 1:
            return r[0]
        else:
            return None

    def result(self, wuid, clientid, files, errorcode=None,
               failedcommand=None, commit=True):
        cursor = self.get_cursor()
        if not self.conn.in_transaction:
            cursor.begin(EXCLUSIVE)
        data = self.get_by_wuid(cursor, wuid)
        if data is None:
            self.commit(commit)
            cursor.close()
            return False
        try:
            self._checkstatus(data, WuStatus.ASSIGNED)
        except StatusUpdateError:
            self.commit(commit)
            cursor.close()
            if data["status"] == WuStatus.CANCELLED:
                global PRINTED_CANCELLED_WARNING
                if not PRINTED_CANCELLED_WARNING:
                    logger.warning("If workunits get cancelled due to timeout "
                            "even though the clients are still processing them, "
                            "consider increasing the tasks.wutimeout parameter or "
                            "decreasing the range covered in each workunit, "
                            "i.e., the tasks.polyselect.adrange or "
                            "tasks.sieve.qrange parameters.")
                    PRINTED_CANCELLED_WARNING = True
            raise
        if DEBUG > 0:
            self.check(data)
        d = {"resultclient": clientid,
             "errorcode": errorcode,
             "failedcommand": failedcommand,
             "timeresult": str(datetime.utcnow())}
        if errorcode is None or errorcode == 0:
           d["status"] = WuStatus.RECEIVED_OK
        else:
            d["status"] = WuStatus.RECEIVED_ERROR
        pk = self.mapper.getpk()
        self._add_files(cursor, files, rowid = data[pk])
        self.mapper.table.update(cursor, d, eq={pk:data[pk]})
        self.commit(commit)
        cursor.close()
        return True

    def verification(self, wuid, ok, commit=True):
        cursor = self.get_cursor()
        if not self.conn.in_transaction:
            cursor.begin(EXCLUSIVE)
        data = self.get_by_wuid(cursor, wuid)
        if data is None:
            self.commit(commit)
            cursor.close()
            return False
        # FIXME: should we do the update by wuid and skip these checks?
        try:
            self._checkstatus(data, [WuStatus.RECEIVED_OK, WuStatus.RECEIVED_ERROR])
        except StatusUpdateError:
            self.commit(commit)
            cursor.close()
            raise
        if DEBUG > 0:
            self.check(data)
        d = {"timeverified": str(datetime.utcnow())}
        d["status"] = WuStatus.VERIFIED_OK if ok else WuStatus.VERIFIED_ERROR
        pk = self.mapper.getpk()
        self.mapper.table.update(cursor, d, eq={pk:data[pk]})
        self.commit(commit)

        cursor.close()
        return True

    def cancel(self, wuid, commit=True):
        self.cancel_by_condition(eq={"wuid": wuid}, commit=commit)

    def cancel_all_available(self, commit=True):
        self.cancel_by_condition(eq={"status": WuStatus.AVAILABLE}, commit=commit)

    def cancel_all_assigned(self, commit=True):
        self.cancel_by_condition(eq={"status": WuStatus.ASSIGNED}, commit=commit)

    def cancel_by_condition(self, commit=True, **conditions):
        self.set_status(WuStatus.CANCELLED, commit=commit, **conditions)

    def set_status(self, status, commit=True, **conditions):
        cursor = self.get_cursor()
        if not self.conn.in_transaction:
            cursor.begin(EXCLUSIVE)
        self.mapper.table.update(cursor, {"status": status}, **conditions)
        self.commit(commit)
        cursor.close()

    def query(self, limit=None, **conditions):
        cursor = self.get_cursor()
        r = self.mapper.where(cursor, limit=limit, **conditions)
        cursor.close()
        return r

    def count(self, **cond):
        cursor = self.get_cursor()
        count = self.mapper.count(cursor, **cond)
        cursor.close()
        return count

    def count_available(self):
        return self.count(eq={"status": WuStatus.AVAILABLE})

    def get_one_result(self):
        r = self.query(limit = 1, eq={"status": WuStatus.RECEIVED_OK})
        if not r:
            r = self.query(limit = 1, eq={"status": WuStatus.RECEIVED_ERROR})
        if not r:
            return None
        else:
            return r[0]
#}

class WuResultMessage(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def get_wu_id(self):
        pass
    @abc.abstractmethod
    def get_output_files(self):
        pass
    @abc.abstractmethod
    def get_stdout(self, command_nr):
        pass
    @abc.abstractmethod
    def get_stdoutfile(self, command_nr):
        pass
    @abc.abstractmethod
    def get_stderr(self, command_nr):
        pass
    @abc.abstractmethod
    def get_stderrfile(self, command_nr):
        pass
    @abc.abstractmethod
    def get_exitcode(self, command_nr):
        pass
    @abc.abstractmethod
    def get_command_line(self, command_nr):
        pass
    @abc.abstractmethod
    def get_host(self):
        pass
    def _read(self, filename, data):
        if not filename is None:
            with open(filename, "rb") as inputfile:
                data = inputfile.read()
        return bytes() if data is None else data
    def read_stdout(self, command_nr):
        """ Returns the contents of stdout of command_nr as a byte string.

        If no stdout was captured, returns the empty byte string.
        """
        return self._read(self.get_stdoutfile(command_nr),
                          self.get_stdout(command_nr))
    def read_stderr(self, command_nr):
        """ Like read_stdout() but for stderr """
        return self._read(self.get_stderrfile(command_nr),
                          self.get_stderr(command_nr))


class ResultInfo(WuResultMessage):
    def __init__(self, record):
        # record looks like this:
        # {'status': 0, 'errorcode': None, 'timeresult': None, 'wuid': 'testrun_polyselect_0-5000',
        #  'wurowid': 1, 'timecreated': '2013-05-23 22:28:08.333310', 'timeverified': None,
        #  'failedcommand': None, 'priority': None, 'wu': "WORKUNIT [..rest of workunit text...] \n",
        #  'assignedclient': None, 'retryof': None, 'timeassigned': None, 'resultclient': None,
        #  'files': None}
        self.record = record

    def __str__(self):
        return str(self.record)
    def get_wu_id(self):
        return self.record["wuid"]

    def get_output_files(self):
        """ Returns the list of output files of this workunit

        Only files that were specified in RESULT lines appear here;
        automatically captured stdout and stderr does not.
        """
        if self.record["files"] is None:
            return []
        files = []
        for f in self.record["files"]:
            if f["type"] == "RESULT":
                files.append(f["path"])
        return files

    def _get_stdio(self, filetype, command_nr):
        """ Get the file location of the stdout or stderr file of the
        command_nr-th command. Used internally.
        """
        if self.record["files"] is None:
            return None
        for f in self.record["files"]:
            if f["type"] == filetype and int(f["command"]) == command_nr:
                return f["path"]
        return None

    def get_stdout(self, command_nr):
        # stdout is always captured into a file, not made available directly
        return None

    def get_stdoutfile(self, command_nr):
        """ Return the path to the file that captured stdout of the
        command_nr-th COMMAND in the workunit, or None if there was no stdout
        output. Note that explicitly redirected stdout that was uploaded via
        RESULT does not appear here, but in get_files()
        """
        return self._get_stdio("stdout", command_nr)

    def get_stderr(self, command_nr):
        # stderr is always captured into a file, not made available directly
        return None

    def get_stderrfile(self, command_nr):
        """ Like get_stdoutfile(), but for stderr """
        return self._get_stdio("stderr", command_nr)

    def get_exitcode(self, command_nr):
        """ Return the exit code of the command_nr-th command """
        if not self.record["failedcommand"] is None \
                and command_nr == int(self.record["failedcommand"]):
            return int(self.record["errorcode"])
        else:
            return 0

    def get_command_line(self, command_nr):
        return None

    def get_host(self):
        return self.record["resultclient"]


class DbListener(patterns.Observable):
    """ Class that queries the Workunit database for available results
    and sends them to its Observers.

    The query is triggered by receiving a SIGUSR1 (the instance subscribes to
    the signal handler relay), or by calling send_result().
    """
    # FIXME: SIGUSR1 handler is not implemented
    def __init__(self, *args, db, **kwargs):
        super().__init__(*args, **kwargs)
        self.wuar = WuAccess(db)

    def send_result(self):
        # Check for results
        r = self.wuar.get_one_result()
        if not r:
            return False
        message = ResultInfo(r)
        was_received = self.notifyObservers(message)
        if not was_received:
            logger.error("Result for workunit %s was not processed by any task. "
                         "Setting it to status CANCELLED", message.get_wu_id())
            self.wuar.cancel(message.get_wu_id())
        return was_received

class IdMap(object):
    """ Identity map. Ensures that DB-backed dictionaries of the same table
    name are instantiated only once.

    Problem: we should also require that the DB is identical, but file names
    are not a unique specifier to a file, and we allow connection objects
    instead of DB file name. Not clear how to test for identity, lacking
    support for this from the sqlite3 module API.
    """
    def __init__(self):
        self.db_dicts = {}

    def make_db_dict(self, db, name):
        key = name
        if not key in self.db_dicts:
            self.db_dicts[key] = DictDbAccess(db, name)
        return self.db_dicts[key]

# Singleton instance of IdMap
idmap = IdMap()

class DbAccess(object):
    """ Base class that lets subclasses create DB-backed dictionaries or
    WuAccess instances on a database whose file name is specified in the db
    parameter to __init__.
    Meant to be used as a cooperative class; it strips the db parameter from
    the parameter list and remembers it in a private variable so that it can
    later be used to open DB connections.
    """

    def __init__(self, *args, db, **kwargs):
        super().__init__(*args, **kwargs)
        self.__db = db

    def get_db_connection(self):
        return self.__db.connect()

    def get_db_filename(self):
        return self.__db.path

    def get_db_uri(self):
        return self.__db.uri

    def make_db_dict(self, name, connection=None):
        if connection is None:
            return idmap.make_db_dict(self.__db, name)
        else:
            return idmap.make_db_dict(connection, name)

    def make_wu_access(self, connection=None):
        if connection is None:
            return WuAccess(self.__db)
        else:
            return WuAccess(connection)

    def make_db_listener(self, connection=None):
        if connection is None:
            return DbListener(db=self.__db)
        else:
            return DbListener(db=connection)


class HasDbConnection(DbAccess):
    """ Gives sub-classes a db_connection attribute which is a database
    connection instance.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.db_connection = self.get_db_connection()


class UsesWorkunitDb(HasDbConnection):
    """ Gives sub-classes a wuar attribute which is WuAccess instance, using
    the sub-classes' shared database connection.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.wuar = self.make_wu_access(self.db_connection)


class DbWorker(DbAccess, threading.Thread):
    """Thread executing WuAccess requests from a given tasks queue"""
    def __init__(self, taskqueue, *args, daemon=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.taskqueue = taskqueue
        if not daemon is None:
            self.daemon = daemon
        self.start()

    def run(self):
        # One DB connection per thread. Created inside the new thread to make
        # sqlite happy
        wuar = self.make_wu_access()
        while True:
            # We expect a 4-tuple in the task queue. The elements of the tuple:
            # a 2-array, where element [0] receives the result of the DB call,
            #  and [1] is an Event variable to notify the caller when the
            #  result is available
            # fn_name, the name (as a string) of the WuAccess method to call
            # args, a tuple of positional arguments
            # kargs, a dictionary of keyword arguments
            (result_tuple, fn_name, args, kargs) = self.taskqueue.get()
            if fn_name == "terminate":
                break
            ev = result_tuple[1]
            # Assign to tuple in-place, so result is visible to caller.
            # No slice etc. here which would create a copy of the array
            try: result_tuple[0] = getattr(wuar, fn_name)(*args, **kargs)
            except Exception as e:
                traceback.print_exc()
            ev.set()
            self.taskqueue.task_done()

class DbRequest(object):
    """ Class that represents a request to a given WuAccess function.
        Used mostly so that DbThreadPool's __getattr__ can return a callable
        that knows which of WuAccess's methods should be called by the
        worker thread """
    def __init__(self, taskqueue, func):
        self.taskqueue = taskqueue
        self.func = func

    def do_task(self, *args, **kargs):
        """Add a task to the queue, wait for its completion, and return the result"""
        ev = threading.Event()
        result = [None, ev]
        self.taskqueue.put((result, self.func, args, kargs))
        ev.wait()
        return result[0]

class DbThreadPool(object):
    """Pool of threads consuming tasks from a queue"""
    def __init__(self, dburi, num_threads=1):
        self.taskqueue = Queue(num_threads)
        self.pool = []
        for _ in range(num_threads):
            worker = DbWorker(self.taskqueue, daemon=True, db=dburi)
            self.pool.append(worker)

    def terminate(self):
        for t in self.pool:
            self.taskqueue.put((None, "terminate", None, None))
        self.wait_completion

    def wait_completion(self):
        """Wait for completion of all the tasks in the queue"""
        self.taskqueue.join()

    def __getattr__(self, name):
        """ Delegate calls to methods of WuAccess to a worker thread.
            If the called method exists in WuAccess, creates a new
            DbRequest instance that remembers the name of the method that we
            tried to call, and returns the DbRequest instance's do_task
            method which will process the method call via the thread pool.
            We need to go through a new object's method since we cannot make
            the caller pass the name of the method to call to the thread pool
            otherwise """
        if hasattr(WuAccess, name):
            task = DbRequest(self.taskqueue, name)
            return task.do_task
        else:
            raise AttributeError(name)


# One entry in the WU DB, including the text with the WU contents
# (FILEs, COMMANDs, etc.) and info about the progress on this WU (when and
# to whom assigned, received, etc.)

    # wuid is the unique wuid of the workunit
    # status is a status code as defined in WuStatus
    # data is the str containing the text of the workunit
    # timecreated is the string containing the date and time of when the WU was added to the db
    # timeassigned is the ... of when the WU was assigned to a client
    # assignedclient is the clientid of the client to which the WU was assigned
    # timeresult is the ... of when a result for this WU was received
    # resultclient is the clientid of the client that uploaded a result for this WU
    # errorcode is the exit status code of the first failed command, or 0 if none failed
    # timeverified is the ... of when the result was marked as verified


if __name__ == '__main__': # {
    import argparse

    queries = {"avail" : ("Available workunits", {"eq": {"status": WuStatus.AVAILABLE}}),
               "assigned": ("Assigned workunits", {"eq": {"status": WuStatus.ASSIGNED}}),
               "receivedok": ("Received ok workunits", {"eq":{"status": WuStatus.RECEIVED_OK}}),
               "receivederr": ("Received with error workunits", {"eq": {"status": WuStatus.RECEIVED_ERROR}}),
               "verifiedok": ("Verified ok workunits", {"eq": {"status": WuStatus.VERIFIED_OK}}),
               "verifiederr": ("Verified with error workunits", {"eq": {"status": WuStatus.VERIFIED_ERROR}}),
               "cancelled": ("Cancelled workunits", {"eq": {"status": WuStatus.CANCELLED}}),
               "all": ("All existing workunits", {})
              }

    use_pool = False

    parser = argparse.ArgumentParser()
    parser.add_argument('-dbfile', help='Name of the database file')
    parser.add_argument('-create', action="store_true",
                        help='Create the database tables if they do not exist')
    parser.add_argument('-add', action="store_true",
                        help='Add new workunits. Contents of WU(s) are '
                        'read from stdin, separated by blank line')
    parser.add_argument('-assign', nargs = 1, metavar = 'clientid',
                        help = 'Assign an available WU to clientid')
    parser.add_argument('-cancel', action="store_true",
                        help = 'Cancel selected WUs')
    parser.add_argument('-expire', action="store_true",
                        help = 'Expire selected WUs')
    # parser.add_argument('-setstatus', metavar = 'STATUS',
    #                    help = 'Forcibly set selected workunits to status (integer)')
    parser.add_argument('-prio', metavar = 'N',
                        help = 'If used with -add, newly added WUs '
                        'receive priority N')
    parser.add_argument('-limit', metavar = 'N',
                        help = 'Limit number of records in queries',
                        default = None)
    parser.add_argument('-result', nargs = 6,
                        metavar = ('wuid', 'clientid', 'filename', 'filepath',
                                   'filetype', 'command'),
                        help = 'Return a result for wu from client')
    parser.add_argument('-test', action="store_true",
                        help='Run some self tests')
    parser.add_argument('-debug', help='Set debugging level')
    parser.add_argument('-setdict', nargs = 4,
                        metavar = ("dictname", "keyname", "type", "keyvalue"),
                        help='Set an entry of a DB-backed dictionary')

    parser.add_argument('-wuid', help="Select workunit with given id",
                        metavar="WUID")
    for arg in queries:
        parser.add_argument('-' + arg, action="store_true", required=False,
                            help="Select %s" % queries[arg][0].lower())
    parser.add_argument('-dump', nargs='?', default = None, const = "all",
                        metavar = "FIELD",
                        help='Dump WU contents, optionally a single field')
    parser.add_argument('-sort', metavar = "FIELD",
                        help='With -dump, sort output by FIELD')
    # Parse command line, store as dictionary
    args = vars(parser.parse_args())
    # print(args)

    dbname = "wudb"
    if args["dbfile"]:
        dbname = args["dbfile"]

    if args["test"]:
        import doctest
        doctest.testmod()

    if args["debug"]:
        DEBUG = int(args["debug"])
    prio = 0
    if args["prio"]:
        prio = int(args["prio"][0])
    limit = args["limit"]

    db = DBFactory('db:sqlite3://%s' % dbname)

    if use_pool:
        db_pool = DbThreadPool(db)
    else:
        db_pool = WuAccess(db)

    if args["create"]:
        db_pool.create_tables()
    if args["add"]:
        s = ""
        wus = []
        for line in sys.stdin:
            if line == "\n":
                wus.append(s)
                s = ""
            else:
                s += line
        if s != "":
            wus.append(s)
        db_pool.create(wus, priority=prio)

    # Functions for queries
    queries_list = []
    for (arg, (msg, condition)) in queries.items():
        if args[arg]:
            queries_list.append([msg, condition])
    if args["wuid"]:
        for wuid in args["wuid"].split(","):
            msg = "Workunit %s" % wuid
            condition = {"eq": {"wuid": wuid}}
            queries_list.append([msg, condition])

    for (msg, condition) in queries_list:
        print("%s: " % msg)
        if not args["dump"]:
            count = db_pool.count(limit=args["limit"], **condition)
            print (count)
        else:
            wus = db_pool.query(limit=args["limit"], **condition)
            if wus is None:
                print("0")
            else:
                print (len(wus))
                if args["sort"]:
                    wus.sort(key=lambda wu: str(wu[args["sort"]]))
                if args["dump"] == "all":
                    print(WuAccess.to_str(wus))
                else:
                    for wu in wus:
                        print(wu[args["dump"]])
        if args["cancel"]:
            print("Cancelling selected workunits")
            db_pool.cancel_by_condition(**condition)
        if args["expire"]:
            print("Expiring selected workunits")
            db_pool.set_status(WuStatus.NEED_RESUBMIT, commit=True, **condition)
        # if args["setstatus"]:
        #    db_pool.set_status(int(args["setstatus"]), **condition)

    # Dict manipulation
    if args["setdict"]:
        (name, keyname, itemtype, keyvalue) = args["setdict"]
        # Type-cast value to the specified type
        value =  getattr(__builtins__, itemtype)(keyvalue)
        dbdict = DictDbAccess(dbname, name)
        dbdict[keyname] = value
        del(dbdict)

    # Functions for testing
    if args["assign"]:
        clientid = args["assign"][0]
        wus = db_pool.assign(clientid)

    if args["result"]:
        result = args["result"]
        db_pool.result(result.wuid, result.clientid, result[2:])

    if use_pool:
        db_pool.terminate()
# }

# Local Variables:
# version-control: t
# End:
