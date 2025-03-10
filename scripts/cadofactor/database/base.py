import abc
import logging
import re
import sys
import threading
import time
import traceback

from cadofactor.database.misc import join3, dict_join3

READONLY = object()
DEFERRED = object()
IMMEDIATE = object()
EXCLUSIVE = object()

logger = logging.getLogger("Database")
logger.setLevel(logging.NOTSET)


class TransactionAborted(RuntimeError):
    def __init__(self, command, values):
        super().__init__(f"{command}")


class pending_transaction(object):
    def __init__(self, *args):
        self._set(*args)

    def _set(self, connection=None, transaction=None):
        self.connection = connection
        self.transaction = transaction
        if connection is not None:
            assert transaction is not None
            self.thread = threading.current_thread()
            self.traceback = traceback.extract_stack()
            assert self.traceback is not None
        else:
            assert transaction is None
            self.thread = None
            self.traceback = None

    def set(self, connection, transaction):
        assert connection is not None
        assert transaction is not None
        if self.connection is None:
            self._set(connection, transaction)
        else:
            incoming = pending_transaction(connection, transaction)
            if self.traceback is not None:
                old_tb_str = "".join(traceback.format_list(self.traceback))
            else:
                old_tb_str = "<no stacktrace???>"
            new_tb_str = "".join(traceback.format_list(incoming.traceback))
            logger.transaction("EXCLUSIVE conflict on %s\n%s",
                               self, old_tb_str)
            logger.transaction("Incoming transaction: %s\n%s",
                               incoming, new_tb_str)

    def release(self):
        self._set()

    def __str__(self):
        if self.connection is None:
            return "<empty>"
        return ", ".join([f"connection 0x{id(self.connection):x}",
                          f"thread 0x{self.thread.ident:x}",
                          f"transaction 0x{id(self.transaction):x}"])

    def __repr__(self):
        return str(self)


class transaction_logbook(dict):
    """
    For each database (and only once per database file name), we want to
    record the connection id, thread id, and transaction id of any
    pending exclusive transaction, and verify that there is only one!
    """
    def __init__(self):
        pass

    def new_db(self, name):
        if name not in self:
            logger.transaction(f"Logging pending transactions for {name}")
            self[name] = pending_transaction()
        return self[name]


pending_transactions = transaction_logbook()


class ConnectionWrapperBase(object, metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def cursor(self):
        pass

    def transaction(self, mode=None):
        """
        harness_transaction(mode, <some lambda>) is preferred.
        """
        return TransactionWrapper(self.cursor(), mode=mode)

    def harness_transaction(self, mode, transaction,
                            *args, replay=True, **kwargs):
        i = 0
        while True:
            try:
                with self.transaction(mode) as cursor:
                    return transaction(cursor, *args, **kwargs)
            except TransactionAborted as e:
                if replay:
                    i += 1
                    logger.warning("rolling back and retrying transaction"
                                   f" (#{i})"
                                   ": thread 0x%08x connection 0x%08x"
                                   " transaction lambda %s",
                                   threading.current_thread().ident,
                                   id(self),
                                   transaction)
                    time.sleep(min(1, 0.125 * 2**i))
                else:
                    raise e


# I wish I knew how to make that inherit from a template argument (which
# would be sqlite3.Cursor or mysql.Cursor). I'm having difficulties to
# grok that syntax though, so let's stay simple and stupid. We'll have a
# *member* which is the cursor object, and so be it.
class CursorWrapperBase(object, metaclass=abc.ABCMeta):
    """
    This class represents a DB cursor and provides convenience functions
    around SQL queries. In particular it is meant to provide an
    (1) an interface to SQL functionality via method calls with parameters,
    and
    (2) hiding some particularities of the SQL variant of the underlying
        DBMS as far as possible
    """

    # This is used in where queries; it converts from named arguments such as
    # "eq" to a binary operator such as "="
    name_to_operator = {"lt": "<",
                        "le": "<=",
                        "eq": "=",
                        "ge": ">=",
                        "gt": ">",
                        "ne": "!=",
                        "like": "like"}

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
        if type(x) is tuple:
            return tuple([self.translations(u) for u in x])
        elif type(x) is list:
            return [self.translations(u) for u in x]
        else:
            v = x
            for a, b in self._string_translations:
                v, nrepl = re.subn(a, b, v)
            return v

    def upgrade_command(self, command):
        # only if needed. Might be useful for selects within
        # transactions.
        return command

    # override in the derived cursor class if needed
    @property
    def parameter_auto_increment(self):
        return "?"

    def __init__(self):
        # This boolean is set to true from within _begin_transaction if
        # we acquire a lock
        self.in_exclusive_transaction = None

    def in_transaction(self):
        return self.connection.in_transaction

    @staticmethod
    def _without_None(d):
        """ Return a copy of the dictionary d, but without entries whose values
            are None """
        return {k[0]: k[1] for k in d.items() if k[1] is not None}

    @staticmethod
    def as_string(d):
        if d is None:
            return ""
        else:
            return ", " + dict_join3(d, sep=", ", op=" AS ")

    def _where_str(self, name, **args):
        where = ""
        values = []
        qm = self.parameter_auto_increment
        for opname in args:
            if args[opname] is None:
                continue
            if where == "":
                where = f" {name} "
            else:
                where += " AND "
            where += join3(args[opname].keys(),
                           post=f" {self.name_to_operator[opname]} {qm}",
                           sep=" AND ")
            values += list(args[opname].values())
        return (where, values)

    def _exec(self, command, values=None):
        """
        Wrapper around self.execute() that prints arguments
        for debugging and retries in case of "database locked" exception
        """

        # FIXME: should be the caller's class name, as _exec could be
        # called from outside this class
        classname = self.__class__.__name__
        parent = sys._getframe(1).f_code.co_name
        command = self.translations(command)
        if self.in_exclusive_transaction is not None:
            command = self.upgrade_command(command)

        command_str = command.replace("?", "%r")

        if values is not None:
            command_str = command_str % tuple(values)

        if self.in_exclusive_transaction is not None:
            logger.transaction("%s.%s():"
                               " connection = 0x%x,"
                               " cursor = 0x%x,"
                               " transaction = 0x%x,"
                               " command = %s",
                               classname, parent,
                               id(self.connection),
                               id(self),
                               id(self.in_exclusive_transaction),
                               command_str)
        else:
            logger.transaction("%s.%s():"
                               " connection = 0x%x,"
                               " cursor = 0x%x,"
                               " command = %s",
                               classname, parent,
                               id(self.connection),
                               id(self),
                               command_str)

        self.try_catch_execute(command, values)

        logger.transaction("%s.%s(): connection = 0x%x, command finished",
                           classname, parent, id(self.connection))

    def try_catch_execute(self, command, values):
        """
        This is a priori overriden by the subclasses
        """
        if values is None:
            self.cursor.execute(command)
        else:
            self.cursor.execute(command, values)

    def _begin_transaction(self, mode=None, transaction=None):
        if mode is None:
            self._exec("BEGIN")
        elif mode is DEFERRED:
            self._exec("BEGIN DEFERRED")
        elif mode is IMMEDIATE:
            self._exec("BEGIN IMMEDIATE")
        elif mode is EXCLUSIVE:
            self._exec("BEGIN EXCLUSIVE")
            self.in_exclusive_transaction = transaction
        else:
            raise TypeError("Invalid mode parameter: %r" % mode)

    def pragma(self, prag):
        self._exec("PRAGMA %s;" % prag)

    def create_table(self, table, layout):
        """
        Creates a table with fields as described in the layout parameter
        """
        command = "CREATE TABLE IF NOT EXISTS %s( %s );" % \
                  (table, ", ".join(" ".join(k) for k in layout))
        self._exec(command)

    def create_index(self, name, table, columns):
        """
        Creates an index with fields as described in the columns list
        """
        # we get so many of these...
        try:
            command = self.translations("CREATE INDEX IF NOT EXISTS")
            command += " %s ON %s( %s );" % (name, table, ", ".join(columns))
            self._exec(command)
        except Exception as e:
            logger.warning(e)
            pass

    def insert(self, table, d):
        """
        Insert a new entry, where d is a dictionary containing the
        field:value pairs. Returns the row id of the newly created entry
        """
        # INSERT INTO table (field_1, field_2, ..., field_n)
        # 	VALUES (value_1, value_2, ..., value_n)

        # Fields is a copy of d but with entries removed that have value None.
        # This is done primarily to avoid having "id" listed explicitly in the
        # INSERT statement, because the DB fills in a new value automatically
        # if "id" is the primary key. But I guess not listing field:NULL items
        # explicitly in an INSERT is a good thing in general
        fields = self._without_None(d)
        fields_str = ", ".join(fields.keys())

        qm = self.parameter_auto_increment
        # sqlformat = "?, ?, ?, " ... "?"
        sqlformat = ", ".join((qm,) * len(fields))
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
        qm = self.parameter_auto_increment
        setstr = join3(d.keys(), post=" = " + qm, sep=", ")
        (wherestr, wherevalues) = self._where_str("WHERE", **conditions)
        command = "UPDATE %s SET %s %s" % (table, setstr, wherestr)
        values = list(d.values()) + wherevalues
        self._exec(command, values)

    def where_query(self,
                    joinsource,
                    col_alias=None,
                    limit=None, offset=None,
                    order=None,
                    **conditions):
        # Table/Column names cannot be substituted,
        # so include in query directly.
        (WHERE, values) = self._where_str("WHERE", **conditions)
        if order is None:
            ORDER = ""
        else:
            if order[1] not in ("ASC", "DESC"):
                raise Exception
            ORDER = " ORDER BY %s %s" % (order[0], order[1])
        if limit is None:
            LIMIT = ""
        else:
            LIMIT = " LIMIT %s" % int(limit)
            if offset is not None:
                LIMIT += " OFFSET %s" % int(offset)
        AS = self.as_string(col_alias)
        command = "SELECT * %s FROM %s %s %s %s" \
                  % (AS, joinsource, WHERE, ORDER, LIMIT)
        return (command, values)

    def where(self, joinsource, col_alias=None,
              limit=None, offset=None,
              order=None,
              values=[], **conditions):
        """ Get a up to "limit" table rows (limit == 0: no limit) where
            the key:value pairs of the dictionary "conditions" are set to the
            same value in the database table """
        (command, newvalues) = self.where_query(joinsource, col_alias, limit,
                                                order, **conditions)
        self._exec(command + ";", values + newvalues)

    def count(self, joinsource, **conditions):
        """
        Count rows where the key:value pairs of the dictionary
        "conditions" are set to the same value in the database table
        """

        # Table/Column names cannot be substituted, so include in query
        # directly.
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
        desc = [k[0] for k in self.cursor.description]

        return [dict(zip(desc, row)) for row in self.cursor.fetchall()]

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
        return r"(?:db:)?(\w+)://(?:(?:(\w+)(?::(.*))?@)?" \
               r"(?:([\w\.]+)|\[([\d:]+)*\])(?::(\d+))?/)?(.*)$"

    def __init__(self, uri, backend_pattern=None):
        self.uri = uri
        foo = re.match(self.general_pattern, uri)
        if not foo:
            raise ValueError("db URI %s does not match regexp %s"
                             % (uri, self.general_pattern))
        self.hostname = foo.group(4)
        self.host_ipv6 = False
        if self.hostname is None:
            self.hostname = foo.group(5)
            if self.hostname is not None:
                self.host_ipv6 = True
        self.backend = foo.group(1)
        if backend_pattern is not None \
           and not re.match(backend_pattern, self.backend):
            raise ValueError("back-end type %s not supported, expected %s"
                             % (self.backend, backend_pattern))

        D = dict(user=foo.group(2),
                 password=foo.group(3),
                 host=self.hostname,
                 port=foo.group(6))

        if D['host'] is None:
            assert D['port'] is None
            assert D['user'] is None
        if D['user'] is None:
            assert D['password'] is None
        if D['port'] is not None:
            D['port'] = int(D['port'])
        D = {k: v for k, v in D.items() if v is not None}

        self.db_connect_args = dict(**D)
        self.db_name = foo.group(7)
        self.talked = False
        logger.info("Database URI is %s" % self.uri_without_credentials)

    @property
    def uri_without_credentials(self):
        text = "db:%s://" % self.backend
        d = self.db_connect_args
        if "host" in d:
            if "user" in d:
                text += "USERNAME"
                if "password" in d:
                    text += ":PASSWORD"
                text += "@"
            if self.host_ipv6:
                text += "[%s]" % d["host"]
            else:
                text += d["host"]
            if "port" in d:
                text += ":%s" % d["port"]
            text += "/"
        text += self.db_name
        return text

    def advertise_connection(self):
        if not self.talked:
            logger.info("Opened connection to database %s" % self.db_name)
            self.talked = True


class TransactionWrapper(object):
    def __init__(self, cursor, mode=None):
        self.cursor = cursor
        self.mode = mode

    def __enter__(self):
        thr_id = threading.current_thread().ident
        conn = self.cursor.connection
        if self.mode == READONLY:
            what = "read-only lookup"
        else:
            what = "transaction"
        logger.transaction(f"thread 0x{thr_id:x}"
                           f" opens {what} 0x{id(self):x}"
                           f" on connection 0x{id(conn):x}"
                           f" [in_transaction: {conn.in_transaction}]")
        if self.mode == READONLY:
            return self.cursor
        if self.mode == EXCLUSIVE:
            conn.pending.set(conn, self)
        self.cursor._begin_transaction(self.mode, self)
        return self.cursor

    def __exit__(self, e_type, e_value, traceback):
        thr_id = threading.current_thread().ident
        conn = self.cursor.connection
        if self.mode == READONLY:
            what = "read-only lookup"
        else:
            what = "transaction"
        if e_type is None:
            conn.commit()
        else:
            logger.transaction(f"thread 0x{thr_id:x}"
                               f" rolls back {what} 0x{id(self):x}"
                               f" on connection 0x{id(conn):x}")
            conn.rollback()
        conn.pending.release()
        self.cursor.close()
        logger.transaction(f"thread 0x{thr_id:x}"
                           f" closed {what} 0x{id(self):x}"
                           f" on connection 0x{id(conn):x}"
                           f" [in_transaction: {conn.in_transaction}]")


def conn_close(conn):
    try:
        logger.transaction("Closing connection %d", id(conn))
        if conn.in_transaction():
            logger.warning("Connection %d being closed while in transaction",
                           id(conn))
        conn.close()
    except Exception:
        pass
