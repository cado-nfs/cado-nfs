#!/usr/bin/env python3

# Make Python 2.7 use the print() syntax from Python 3
from __future__ import print_function

import sys
import sqlite3
import threading
import traceback
from datetime import datetime
from workunit import Workunit
if sys.version_info.major == 3:
    from queue import Queue
else:
    from Queue import Queue

debug = 1

def diag(level, text, var = None):
    if debug > level:
        if var is None:
            print (text, file=sys.stderr)
        else:
            print (text + str(var), file=sys.stderr)

def join3(l, op = None, pre = None, post = None, sep = ", "):
    """ For a list l = ('a', 'b', 'c'), string pre = "+", 
    string post = "-", and set = ", ", 
    returns the string '+a-, +b-, +c-' 
    If any parameter is None, it is interpreted as the empty string """
    if pre is None:
        pre = ""
    if post is None:
        post = ""
    if sep is None:
        sep = "";
    if op is None:
        op = ""
    return sep.join([pre + op.join(list(k)) + post for k in l])

def dict_join3(d, sep=None, op=None, pre=None, post=None):
    """ For a dictionary {"a": "1", "b": "2"}, sep = "," op = "=",
    pre="-" and post="+", returns "-a=1+,-b=2+"
    If any parameter is None, it is interpreted as the empty string """
    return join3(d.items(), sep=sep, op=op, pre=pre, post=post)

# Dummy class for defining "constants"
class WuStatus:
    AVAILABLE = 0
    ASSIGNED = 1
    RECEIVED_OK = 2
    RECEIVED_ERROR = 3
    VERIFIED_OK = 4
    VERIFIED_ERROR = 5
    CANCELLED = 6

    @classmethod
    def check(cls, status):
        """ Check whether status is equal to one of the constants """
        assert status in (WuStatus.AVAILABLE, WuStatus.ASSIGNED, WuStatus.RECEIVED_OK, 
            WuStatus.RECEIVED_ERROR, WuStatus.VERIFIED_OK, WuStatus.VERIFIED_ERROR, 
            WuStatus.CANCELLED)

# If we try to update the status in any way other than progressive 
# (AVAILABLE -> ASSIGNED -> ...), we raise this exception
class StatusUpdateError(Exception):
    pass

class MyCursor(sqlite3.Cursor):
    """ This class represents a DB cursor and provides convenience functions 
        around SQL queries. In particular it is meant to provide an  
        (1) an interface to SQL functionality via method calls with parameters, 
        and 
        (2) hiding some particularities of the SQL variant of the underlying 
            DBMS as far as possible """
        
    # This is used in where queries; it converts from named arguments such as 
    # "eq" to a binary operator such as "="
    name_to_operator = {"lt": "<", "le": "<=", "eq": "=", "ge": ">=", "gt" : ">", "ne": "!="}
    
    def __init__(self, conn):
        # Enable foreign key support
        super(MyCursor, self).__init__(conn)
        self._exec("PRAGMA foreign_keys = ON;")

    @staticmethod
    def _fieldlist(l, r = "=", s = ", "):
        """ For a list l = ('a', 'b', 'c') returns the string 'a = ?, b = ?, c = ?',
            or with a different string r in place of the "=", or a different 
            string s in place of the ", " """
        return s.join([k + " " + r + " ?" for k in l])

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
    
    @classmethod
    def where_str(cls, name, **args):
        where = ""
        values = []
        for opname in args:
            if args[opname] is None:
                continue
            if where == "":
                where = " " + name + " "
            else:
                where = where + " AND "
            where = where + cls._fieldlist(args[opname].keys(), cls.name_to_operator[opname], s = " AND ")
            values = values + list(args[opname].values())
        return (where, values)

    def _exec(self, command, values = None):
        """ Wrapper around self.execute() that prints arguments 
            for debugging and retries in case of "database locked" exception """
        if debug > 1:
            # FIXME: should be the caller's class name, as _exec could be 
            # called from outside this class
            classname = self.__class__.__name__
            parent = sys._getframe(1).f_code.co_name
            diag (1, classname + "." + parent + "(): command = " + command)
            if not values is None:
                diag (1, classname + "." + parent + "(): values = ", values)
        i = 0
        while True:
            try:
                if values is None:
                    self.execute(command)
                else:
                    self.execute(command, values)
                break
            except sqlite3.OperationalError as e:
                if i == 10 or str(e) != "database is locked":
                    raise

    def create_table(self, table, layout):
        """ Creates a table with fields as described in the layout parameter """
        command = "CREATE TABLE IF NOT EXISTS " + table + \
            "( " + join3(layout, op=" ", sep=", ") + " );"
        self._exec (command)
    
    def create_index(self, table, d):
        """ Creates an index with fields as described in the d dictionary """
        for (name, columns) in d.items():
            column_names = [col[0] for col in columns]
            command = "CREATE INDEX IF NOT EXISTS " + name + " ON " + \
                table + "( " + ", ".join(column_names) + " );"
            self._exec (command)
    
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

        sqlformat = ", ".join(("?",) * len(fields)) # sqlformat = "?, ?, ?, " ... "?"
        command = "INSERT INTO " + table + \
            " (" + ", ".join(fields.keys()) + ") VALUES (" + sqlformat + ");"
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
        setstr = " SET " + self.__class__._fieldlist(d.keys())
        setvalues = d.values()
        (wherestr, wherevalues) = self.__class__.where_str("WHERE", **conditions)
        command = "UPDATE " + table + setstr + wherestr
        values = list(setvalues) + wherevalues
        self._exec(command, values)
    
    def where(self, joinsource, col_alias = None, limit = None, order = None, 
              **conditions):
        """ Get a up to "limit" table rows (limit == 0: no limit) where 
            the key:value pairs of the dictionary d are set to the same 
            value in the database table """

        # Table/Column names cannot be substituted, so include in query directly.
        (WHERE, values) = self.__class__.where_str("WHERE", **conditions)

        if order is None:
            ORDER = ""
        else:
            if not order[1] in ("ASC", "DESC"):
                raise Exception
            ORDER = " ORDER BY " + str(order[0]) + " " + order[1]

        if limit is None:
            LIMIT = ""
        else:
            LIMIT = " LIMIT " + str(int(limit))

        AS = self.as_string(col_alias);

        command = "SELECT *" + AS + " FROM " + joinsource + WHERE + \
            ORDER + LIMIT + ";"
        self._exec(command, values)
        

    def where_as_dict(self, joinsource, col_alias = None, limit = None, 
                      order = None, **conditions):
        self.where(joinsource, col_alias=col_alias, limit=limit, 
                      order=order, **conditions)
        # cursor.description is a list of lists, where the first element of 
        # each inner list is the column name
        result = []
        desc = [k[0] for k in self.description]
        row = self.fetchone()
        while row is not None:
            diag (2, "MyCursor.where_as_dict(): row = ", row)
            result.append(dict(zip(desc, row)))
            row = self.fetchone()
        return result


class DbTable(object):
    """ A class template defining access methods to a database table """
    def __init__(self):
        self.tablename = type(self).name
        self.fields = type(self).fields
        self.primarykey = type(self).primarykey
        self.references = type(self).references

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
            fk = (r.primarykey, "INTEGER", "REFERENCES " + 
                  r.getname() + " (" + r.primarykey + ") ")
            fields.append(fk)
        cursor.create_table(self.tablename, fields)
        cursor.create_index(self.tablename, self.index)

    def insert(self, cursor, values, foreign=None):
        """ Insert a new row into this table. The column:value pairs are 
            specified key:value pairs of the dictionary d. 
            The database's row id for the new entry is stored in d[primarykey] """
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

    def where(self, cursor, limit = None, order = None, **conditions):
        assert order is None or order[0] in self._get_colnames()
        return cursor.where_as_dict(self.tablename, limit=limit, 
                                    order=order, **conditions)


class WuTable(DbTable):
    name = "workunits"
    fields = (
        ("wurowid", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"), 
        ("wuid", "TEXT", "UNIQUE NOT NULL"), 
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
        ("retryof", "TEXT", ""),
        ("priority", "INTEGER", "")
    )
    primarykey = fields[0][0]
    references = None
    index = {"wuidindex": (fields[1],), "statusindex" : (fields[2],)}

class FilesTable(DbTable):
    name = "files"
    fields = (
        ("filesrowid", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"), 
        ("filename", "TEXT", ""), 
        ("path", "TEXT", "UNIQUE NOT NULL")
    )
    primarykey = fields[0][0]
    references = WuTable()
    index = {"wuindex": (fields[1],)}


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
    
    def where(self, cursor, limit = None, order = None, **cond):
        pk = self.getpk()
        joinsource = self.table.name
        for s in self.subtables.keys():
            # FIXME: this probably breaks with more than 2 tables
            joinsource = joinsource + " LEFT JOIN " + \
                self.subtables[s].getname() + \
                " USING ( " + pk + " )"
        # FIXME: don't get result rows as dict! Leave as tuple and
        # take them apart positionally
        rows = cursor.where_as_dict(joinsource, limit=limit, order=order, 
                                    **cond)
        wus = []
        for r in rows:
            # Collapse rows with identical primary key
            if len(wus) == 0 or r[pk] != wus[-1][pk]:
                wus.append(self.table.dictextract(r))
                for s in self.subtables.keys():
                    wus[-1][s] = None

            for (sn, sm) in self.subtables.items():
                spk = sm.getpk()
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
        {"id": int, "wuid": string, "filename": string, "path": string
        Operations on instances of WuAcccess are directly carried 
        out on the database persistent storage, i.e., they behave kind 
        of as if the WuAccess instance were itself a persistent 
        storage device """
    
    def __init__(self, conn):
        self.conn = conn
        self.mapper = Mapper(WuTable(), {"files": FilesTable()})

    @staticmethod
    def to_str(wus):
        s = ""
        for wu in wus:
            s = s + "Workunit " + str(wu["wuid"]) + ":\n"
            for (k,v) in wu.items():
                if k != "wuid" and k != "files":
                    s = s + "  " + k + ": " + repr(v) + "\n"
            if "files" in wu:
                s = s + "  Associated files:\n"
                if wu["files"] is None:
                    s = s + "    None\n"
                else:
                    for f in wu["files"]:
                        s = s + "    " + str(f) + "\n"
        return s

    @staticmethod
    def _checkstatus(wu, status):
        diag (2, "WuAccess._checkstatus(" + str(wu) + ", " + str(status) + ")")
        if wu["status"] != status:
            msg = "WU " + str(wu["wuid"]) + " has status " + str(wu["status"]) \
                + ", expected " + str(status)
            diag (0, "WuAccess._checkstatus(): " + msg)
            # FIXME: this raise has no effect other than returning from the 
            # method. Why?
            # raise wudb.StatusUpdateError(msg)
            raise Exception(msg)

    def check(self, data):
        status = data["status"]
        WuStatus.check(status)
        wu = Workunit(data["wu"])
        assert wu.get_id() == data["wuid"]
        if status > WuStatus.RECEIVED_ERROR:
            return
        if status == WuStatus.RECEIVED_ERROR:
            assert data["errorcode"] != 0
            return
        if status == WuStatus.RECEIVED_OK:
            assert data["errorcode"] == 0
            return
        assert data["errorcode"] is None
        assert data["timeresult"] is None
        assert data["resultclient"] is None
        if status == WuStatus.ASSIGNED:
            return
        assert data["timeassigned"] is None
        assert data["assignedclient"] is None
        if status == WuStatus.AVAILABLE:
            return
        assert data["timecreated"] is None
        # etc.
    
    # Here come the application-visible functions that implement the 
    # "business logic": creating a new workunit from the text of a WU file,
    # assigning it to a client, receiving a result for the WU, marking it as
    # verified, or marking it as cancelled

    def add_files(self, cursor, files, wuid = None, rowid = None):
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
        d = ({"filename": f[0], "path": f[1]} for f in files)
        # These two should behave identically
        if True:
            self.mapper.insert(cursor, [{pk:rowid, "files": d},])
        else:
            self.mapper.subtables["files"].insert(cursor, d, foreign=rowid)

    def create_tables(self):
        cursor = self.conn.cursor(MyCursor)
        cursor._exec("PRAGMA journal_mode=WAL;")
        self.mapper.create(cursor)
        self.conn.commit()
        cursor.close()

    def create1(self, cursor, wutext, priority = None):
        d = {
            "wuid": Workunit(wutext).get_id(),
            "wu": wutext,
            "status": WuStatus.AVAILABLE,
            "timecreated": str(datetime.now())
            }
        if not priority is None:
            d["priority"] = priority
        # Insert directly into wu table
        self.mapper.table.insert(cursor, d)

    def create(self, wus, priority = None):
        """ Create new workunits from wus which contains the texts of the 
            workunit files """
        cursor = self.conn.cursor(MyCursor)
        if isinstance(wus, str):
            self.create1(cursor, wus, priority)
        else:
            for wu in wus:
                self.create1(cursor, wu, priority)
        self.conn.commit()
        cursor.close()

    def assign(self, clientid):
        """ Finds an available workunit and assigns it to clientid.
            Returns the text of the workunit, or None if no available 
            workunit exists """
        cursor = self.conn.cursor(MyCursor)
        r = self.mapper.table.where(cursor, limit = 1, 
                                    order=("priority", "DESC"), 
                                    eq={"status": WuStatus.AVAILABLE})
        assert len(r) <= 1
        if len(r) == 1:
            self._checkstatus(r[0], WuStatus.AVAILABLE)
            if debug > 0:
                self.check(r[0])
            d = {"status": WuStatus.ASSIGNED, 
                 "assignedclient": clientid,
                 "timeassigned": str(datetime.now())}
            pk = self.mapper.getpk()
            self.mapper.table.update(cursor, d, eq={pk:r[0][pk]})
            self.conn.commit()
            cursor.close()
            return r[0]["wu"]
        else:
            cursor.close()
            return None

    def get_by_wuid(self, cursor, wuid):
        r = self.mapper.where(cursor, eq={"wuid": wuid})
        assert len(r) <= 1
        if len(r) == 1:
            return r[0]
        else:
            return None

    def result(self, wuid, clientid, files, errorcode = None, 
               failedcommand = None):
        cursor = self.conn.cursor(MyCursor)
        data = self.get_by_wuid(cursor, wuid)
        if data is None:
            cursor.close()
            return False
        self._checkstatus(data, WuStatus.ASSIGNED)
        if debug > 0:
            self.check(data)
        d = {"resultclient": clientid,
             "errorcode": errorcode,
             "failedcommand": failedcommand, 
             "timeresult": str(datetime.now())}
        if errorcode is None or errorcode == 0:
           d["status"] = WuStatus.RECEIVED_OK
        else:
            d["status"] = WuStatus.RECEIVED_ERROR
        pk = self.mapper.getpk()
        self.mapper.table.update(cursor, d, eq={pk:data[pk]})
        self.add_files(cursor, files, rowid = data[pk])
        self.conn.commit()
        cursor.close()

    def verification(self, wuid, ok):
        cursor = self.conn.cursor(MyCursor)
        data = self.get_by_wuid(cursor, wuid)
        if data is None:
            cursor.close()
            return False
        # FIXME: should we do the update by wuid and skip these checks?
        self._checkstatus(data, WuStatus.RECEIVED_OK)
        if debug > 0:
            self.check(data)
        d = {["timeverified"]: str(datetime.now())}
        if ok:
            d["status"] = WuStatus.VERIFIED_OK
        else:
            d["status"] = WuStatus.VERIFIED_ERROR
        pk = self.mapper.getpk()
        self.mapper.table.update(cursor, d, eq={pk:data[pk]})
        self.conn.commit()
        cursor.close()

    def cancel(self, wuid):
        cursor = self.conn.cursor(MyCursor)
        d = {"status": WuStatus.CANCELLED}
        self.mapper.table.update(cursor, d, eq={"wuid": wuid})
        self.conn.commit()
        cursor.close()

    def query(self, limit = None, **conditions):
        cursor = self.conn.cursor(MyCursor)
        r = self.mapper.where(cursor, limit=limit, **conditions)
        cursor.close()
        return r


class DbWorker(threading.Thread):
    """Thread executing WuAccess requests from a given tasks queue"""
    def __init__(self, dbfilename, taskqueue):
        threading.Thread.__init__(self)
        self.dbfilename = dbfilename
        self.taskqueue = taskqueue
        self.start()
    
    def run(self):
        # One DB connection per thread. Created inside the new thread to make
        # sqlite happy
        self.connection = sqlite3.connect(self.dbfilename)
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
            wuar = WuAccess(self.connection)
            # Assign to tuple in-place, so result is visible to caller. 
            # No slice etc. here which would create a copy of the array
            try: result_tuple[0] = getattr(wuar, fn_name)(*args, **kargs)
            except Exception as e: 
                traceback.print_exc()
            ev.set()
            self.taskqueue.task_done()
        self.connection.close()

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
    def __init__(self, dbfilename, num_threads = 1):
        self.taskqueue = Queue(num_threads)
        self.pool = []
        for _ in range(num_threads): 
            self.pool.append(DbWorker(dbfilename, self.taskqueue))

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

    use_pool = False

    parser = argparse.ArgumentParser()
    parser.add_argument('-create', action="store_true", required=False, 
                        help='Create the database tables if they do not exist')
    parser.add_argument('-add', action="store_true", required=False, 
                        help='Add new work units. Contents of WU(s) are ' + 
                        'read from stdin, separated by blank line')
    parser.add_argument('-assign', required = False, nargs = 1, 
                        metavar = 'clientid', 
                        help = 'Assign an available WU to clientid')
    parser.add_argument('-prio', required = False, nargs = 1, metavar = 'N', 
                        help = 'If used with -add, newly added WUs ' + 
                        'receive priority N')
    parser.add_argument('-result', required = False, nargs = 4, 
                        metavar = ('clientid', 'wuid', 'filename', 'filepath'), 
                        help = 'Return a result for wu from client')

    for arg in ("avail", "assigned", "receivedok", "receivederr", "all", 
                "dump"):
        parser.add_argument('-' + arg, action="store_true", required = False)
    for arg in ("dbname", "debug"):
        parser.add_argument('-' + arg, required = False, nargs = 1)
    # Parse command line, store as dictionary
    args = vars(parser.parse_args())
    # print(args)

    dbname = "wudb"
    if args["dbname"]:
        dbname = args["dbname"]

    if args["debug"]:
        debug = int(args["debug"][0])
    prio = 0
    if args["prio"]:
        prio = int(args["prio"][0])
    
    if use_pool:
        db_pool = DbThreadPool(dbname)
    else:
        conn = sqlite3.connect(dbname)
        db_pool = WuAccess(conn)
    
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
                s = s + line
        if s != "":
            wus.append(s)
        db_pool.create(wus, priority=prio)
    # Functions for queries
    if args["avail"]:
        wus = db_pool.query(eq={"status": WuStatus.AVAILABLE})
        print("Available workunits: ")
        if wus is None:
            print(wus)
        else:
            print (len(wus))
            if args["dump"]:
                print(WuAccess.to_str(wus))
    if args["assigned"]:
        wus = db_pool.query(eq={"status": WuStatus.ASSIGNED})
        print("Assigned workunits: ")
        if wus is None:
            print(wus)
        else:
            print (len(wus))
            if args["dump"]:
                print(WuAccess.to_str(wus))
    if args["receivedok"]:
        wus = db_pool.query(eq={"status": WuStatus.RECEIVED_OK})
        print("Received ok workunits: ")
        if wus is None:
            print(wus)
        else:
            print (len(wus))
            if args["dump"]:
                print(WuAccess.to_str(wus))
    if args["receivederr"]:
        wus = db_pool.query(eq={"status": WuStatus.RECEIVED_ERROR})
        print("Received with error workunits: ")
        if wus is None:
            print(wus)
        else:
            print (len(wus))
            if args["dump"]:
                print(WuAccess.to_str(wus))
    if args["all"]:
        wus = db_pool.query()
        print("Existing workunits: ")
        if wus is None:
            print(wus)
        else:
            print (len(wus))
            if args["dump"]:
                print(WuAccess.to_str(wus))
    # Functions for testing
    if args["assign"]:
        clientid = args["assign"][0]
        wus = db_pool.assign(clientid)

    if args["result"]:
        (wuid, clientid, filename, filepath) = args["result"]
        db_pool.result(wuid, clientid, [(filename, filepath)])
    
    if use_pool:
        db_pool.terminate()
    else:
        conn.close()
# }

# Local Variables:
# version-control: t
# End:
