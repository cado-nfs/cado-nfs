#!/usr/bin/env python3

# TODO:
# FILES table: OBSOLETE column
#     OBSOLETE says that this file was replaced by a newer version, for
#     example checkrels may want to create a new file with only part of
#     the output.  Should this be a pointer to the file that replaced the
#     obsolete one? How to signal a file that is obsolete, but not
#     replaced by anything?
#
#     If one file is replaced by several (say, due to a data corruption
#     in the middle), we need a 1:n relationship.
#
#     If several files are replaced by one (merge), we need n:1.
#
#     What do? Do we really want an n:n relationship here? Disallow
#     fragmenting files, or maybe simply not track it in the DB if we do?

# FILES table: CHECKSUM column
#     We need a fast check that the information stored in the DB still
#     accurately reflects the file system contents. The test should also
#     warn about files in upload/ which are not listed in DB


import sys

import collections
import abc
from datetime import datetime
import time

from cadofactor.workunit import Workunit

# we need the cadologger import because it adds the TRANSACTION log
# level. It's ugly
from cadofactor import patterns, cadologger  # noqa: F401

from cadofactor.database import DBFactory, DbTable, DbAccess
from cadofactor.database import READONLY, EXCLUSIVE
from cadofactor.database import DictDbAccess

# XXX we should get rid of these implementation details
from cadofactor.database.sqlite3 import DB_SQLite
from cadofactor.database.base import conn_close

from cadofactor.database.base import logger

DEBUG = 0
PRINTED_CANCELLED_WARNING = False


# Dummy class for defining "constants" with reverse lookup
STATUS_NAMES = ["AVAILABLE", "ASSIGNED", "NEED_RESUBMIT", "RECEIVED_OK",
                "RECEIVED_ERROR", "VERIFIED_OK",
                "VERIFIED_ERROR", "CANCELLED"]
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


class WuTable(DbTable):
    tablename = "workunits"
    fields = (("wurowid", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"),
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
              ("retryof", "INTEGER", "REFERENCES %s(wurowid)" % tablename),
              ("priority", "INTEGER", "")
              )
    primarykey = fields[0][0]
    references = None
    index = {"wuid": (fields[1][0],),
             "submitter": (fields[2][0],),
             "priority": (fields[14][0],),
             "status": (fields[3][0],)
             }


class FilesTable(DbTable):
    tablename = "files"
    fields = (("filesrowid", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"),
              ("filename", "TEXT", ""),
              ("path", "VARCHAR(512)", "UNIQUE NOT NULL"),
              ("type", "TEXT", ""),
              ("command", "INTEGER", "")
              )
    primarykey = fields[0][0]
    references = WuTable()
    index = {}


class Mapper(object):
    """ This class translates between application objects, i.e., Python
        directories, and the relational data layout in an SQL DB, i.e.,
        one or more tables which possibly have foreign key relationships
        that map to hierarchical data structures. For now, only one
        foreign key / subdirectory."""

    def __init__(self, table, subtables=None):
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
            if t is not None:
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
            if pk not in wuc:
                self.table.insert(cursor, wuc, foreign=foreign)
                # Copy primary key into caller's data
                wu[pk] = wuc[pk]
            for subtable_name in sub_dict.keys():
                self.subtables[subtable_name].insert(
                    cursor, sub_dict[subtable_name], foreign=wu[pk])

    def update(self, cursor, wus):
        pk = self.getpk()
        for wu in wus:
            assert wu[pk] is not None
            wuc = wu.copy()
            sub_dict = self.__sub_dict(wuc)
            rowid = wuc.pop(pk, None)
            if rowid:
                # hmm. We had a weird {wp: rowid} here, which doesn't
                # make sense (there's no wp anywhere). This had been
                # around ever since that line was first committed. Maybe
                # it's pk?
                self.table.update(cursor, wuc, {pk: rowid})
            for s in self.subtables.keys():
                self.subtables[s].update(cursor, sub_dict[s])

    def count(self, cursor, **cond):
        joinsource = self.table.tablename
        return cursor.count(joinsource, **cond)

    def where(self, cursor, limit=None, order=None, **cond):
        # We want:
        # SELECT * FROM (SELECT * from workunits
        # WHERE status = 2 LIMIT 1) LEFT JOIN files USING ( wurowid );
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
                if spk in r and r[spk] is not None:
                    if wus[-1][sn] is None:
                        # If this sub-array is empty, init it
                        wus[-1][sn] = [sm.table.dictextract(r)]
                    elif r[spk] != wus[-1][sn][-1][spk]:
                        # If not empty, and primary key of sub-table is not
                        # same as in previous entry, add it
                        wus[-1][sn].append(sm.table.dictextract(r))
        return wus


class WuAccess(object):
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

        def t(cursor):
            if isinstance(cursor, DB_SQLite.CursorWrapper):
                cursor.pragma("foreign_keys = ON")

        self.conn.harness_transaction(READONLY, t)
        self.mapper = Mapper(WuTable(), {"files": FilesTable()})

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
            for k, v in wu.items():
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
        # logger.debug("WuAccess._checkstatus(%s, %s)", wu, status)
        wu_status = wu["status"]
        if isinstance(status, collections.abc.Container):
            ok = wu_status in status
        else:
            ok = wu_status == status
        if ok:
            return
        msg = "Workunit %s has status %s (%s), expected %s (%s)" \
              % (wu["wuid"], wu_status, WuStatus.get_name(wu_status),
                 status, WuStatus.get_name(status))

        if status is WuStatus.ASSIGNED and wu_status is WuStatus.CANCELLED:
            logger.warning("WuAccess._checkstatus(): %s,"
                           " presumably timed out", msg)
            raise StatusUpdateError(msg)
        elif status is WuStatus.ASSIGNED \
                and wu_status is WuStatus.NEED_RESUBMIT:  # noqa: E127
            logger.warning("WuAccess._checkstatus(): %s,"
                           " manually expired", msg)
            raise StatusUpdateError(msg)
        else:
            logger.error("WuAccess._checkstatus(): %s", msg)
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
        assert wuid is not None or rowid is not None
        assert wuid is None or rowid is None
        # FIXME: allow selecting row to update directly via wuid, without
        # doing query for rowid first
        pk = self.mapper.getpk()
        if rowid is None:
            wu = self._get_by_wuid(cursor, wuid)
            if wu:
                rowid = wu[pk]
            else:
                return False
        colnames = ("filename", "path", "type", "command")
        # zipped length is that of shortest list, so "command" is optional
        d = (dict(zip(colnames, f)) for f in files)
        # These two should behave identically
        if True:
            self.mapper.insert(cursor, [{pk: rowid, "files": d},])
        else:
            self.mapper.subtables["files"].insert(cursor, d, foreign=rowid)

    def create_tables(self):
        self.conn.harness_transaction(EXCLUSIVE, self.mapper.create)

    def create(self, wus, priority=None, commit=True):
        """ Create new workunits from wus which contains the texts of the
            workunit files """
        if isinstance(wus, Workunit):
            wus = [wus]

        for wu in wus:
            d = {
                "wuid": wu.get_id(),
                "wu": str(wu),
                "status": WuStatus.AVAILABLE,
                "timecreated": str(datetime.utcnow())
                }
            if priority is not None:
                d["priority"] = priority

            # Insert directly into wu table
            self.conn.harness_transaction(EXCLUSIVE,
                                          self.mapper.table.insert,
                                          # cursor is implicitly added
                                          d)

    def assign(self, clientid, commit=True, timeout_hint=None):
        """ Finds an available workunit and assigns it to clientid.
            Returns workunit object, or None if no available
            workunit exists """

        def transaction(cursor):
            r = self.mapper.table.where(cursor, limit=1,
                                        eq={"status": WuStatus.AVAILABLE})

            # This alternative "priority" stuff is the root cause for the
            # server taking time to hand out WUs when the count of
            # available WUs drops to zero.  (introduced in 90ae4beb7 --
            # it's an optional-and-never-used feature anyway)
            # r = self.mapper.table.where(cursor, limit=1,
            #                             order=("priority", "DESC"),
            #                             eq={"status": WuStatus.AVAILABLE})

            if not len(r):
                return None

            self._checkstatus(r[0], WuStatus.AVAILABLE)

            if DEBUG > 0:
                self.check(r[0])

            d = {"status": WuStatus.ASSIGNED,
                 "assignedclient": clientid,
                 "timeassigned": str(datetime.utcnow())
                 }
            pk = self.mapper.getpk()

            self.mapper.table.update(cursor, d, eq={pk: r[0][pk]})
            result = Workunit(r[0]["wu"])
            if timeout_hint:
                result['deadline'] = time.time() + float(timeout_hint)

            return result

        return self.conn.harness_transaction(EXCLUSIVE,
                                             transaction)

    def _get_by_wuid(self, cursor, wuid):
        r = self.mapper.where(cursor, eq={"wuid": wuid})
        assert len(r) <= 1
        if len(r) == 1:
            return r[0]
        else:
            return None

    def result(self, wuid, clientid, files,
               errorcode=None,
               failedcommand=None, commit=True):

        def t(cursor):
            data = self._get_by_wuid(cursor, wuid)
            if data is None:
                return False
            try:
                self._checkstatus(data, WuStatus.ASSIGNED)
            except StatusUpdateError:
                if data["status"] == WuStatus.CANCELLED:
                    global PRINTED_CANCELLED_WARNING
                    if not PRINTED_CANCELLED_WARNING:
                        logger.warning(
                            "If workunits get cancelled due to timeout"
                            " even though the clients are"
                            " still processing them,"
                            " consider increasing the tasks.wutimeout"
                            " parameter or decreasing the range covered"
                            " in each workunit, i.e.,"
                            " the tasks.polyselect.adrange"
                            " or tasks.sieve.qrange parameters.")
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
            self._add_files(cursor, files, rowid=data[pk])
            self.mapper.table.update(cursor, d, eq={pk: data[pk]})
            return True

        return self.conn.harness_transaction(EXCLUSIVE, t)

    def verification(self, wuid, ok, commit=True):

        def t(cursor):
            data = self._get_by_wuid(cursor, wuid)
            if data is None:
                return False
            # FIXME: should we do the update by wuid and skip these checks?
            self._checkstatus(data, [WuStatus.RECEIVED_OK,
                                     WuStatus.RECEIVED_ERROR])
            if DEBUG > 0:
                self.check(data)
            d = {"timeverified": str(datetime.utcnow())}
            if ok:
                d["status"] = WuStatus.VERIFIED_OK
            else:
                d["status"] = WuStatus.VERIFIED_ERROR
            pk = self.mapper.getpk()
            self.mapper.table.update(cursor, d, eq={pk: data[pk]})
            return True

        return self.conn.harness_transaction(EXCLUSIVE, t)

    def cancel(self, wuid, commit=True):
        self.cancel_by_condition(eq={"wuid": wuid},
                                 commit=commit)

    def cancel_all_available(self, commit=True):
        self.cancel_by_condition(eq={"status": WuStatus.AVAILABLE},
                                 commit=commit)

    def cancel_all_assigned(self, commit=True):
        self.cancel_by_condition(eq={"status": WuStatus.ASSIGNED},
                                 commit=commit)

    def cancel_by_condition(self, commit=True, **conditions):
        self.set_status(WuStatus.CANCELLED,
                        commit=commit,
                        **conditions)

    def set_status(self, status, commit=True, **conditions):
        self.conn.harness_transaction(EXCLUSIVE,
                                      self.mapper.table.update,
                                      # cursor is implicitly added
                                      {"status": status},
                                      **conditions)

    def query(self, limit=None, **conditions):
        return self.conn.harness_transaction(READONLY,
                                             self.mapper.where,
                                             # cursor is implicitly added
                                             limit=limit,
                                             **conditions)

    def count(self, **cond):
        return self.conn.harness_transaction(READONLY,
                                             self.mapper.count,
                                             # cursor is implicitly added
                                             **cond)

    def count_available(self):
        return self.count(eq={"status": WuStatus.AVAILABLE})

    def get_one_result(self):
        r = self.query(limit=1, eq={"status": WuStatus.RECEIVED_OK})
        if not r:
            r = self.query(limit=1, eq={"status": WuStatus.RECEIVED_ERROR})
        if r:
            return r[0]


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
        if filename is not None:
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
        # {'status': 0, 'errorcode': None, 'timeresult': None,
        #  'wuid': 'testrun_polyselect_0-5000',
        #  'wurowid': 1, 'timecreated': '2013-05-23 22:28:08.333310',
        #  'timeverified': None,
        #  'failedcommand': None, 'priority': None,
        #  'wu': "WORKUNIT [..rest of workunit text...] \n",
        #  'assignedclient': None, 'retryof': None, 'timeassigned': None,
        #  'resultclient': None,
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
            if f["type"].startswith("RESULT"):
                files.append(f["path"])
        return files

    def _get_stdio(self, filetype, command_nr):
        """ Get the file location of the stdout or stderr file of the
        command_nr-th command. Used internally.
        """
        if self.record["files"] is None:
            return None
        for f in self.record["files"]:
            if f["type"].startswith(filetype) \
               and int(f["command"]) == command_nr:
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
        if self.record["failedcommand"] is not None \
                and command_nr == int(self.record["failedcommand"]):
            return int(self.record["errorcode"])
        else:
            return 0

    def get_command_line(self, command_nr):
        return None

    def get_host(self):
        return self.record["resultclient"]


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
        self.wuar = WuAccess(self.db_connection)


class ListensToWorkunitDb(UsesWorkunitDb, patterns.Observable):
    """ Class that queries the Workunit database for available results
    and sends them to its Observers.

    The query is triggered by receiving a SIGUSR1 (the instance subscribes to
    the signal handler relay), or by calling send_result().
    """
    # FIXME: SIGUSR1 handler is not implemented
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # already in UsesWorkunitDb
        # self.wuar = WuAccess(db)

    def send_result(self):
        # Check for results
        r = self.wuar.get_one_result()
        if not r:
            return False
        message = ResultInfo(r)
        was_received = self.notifyObservers(message)
        if not was_received:
            logger.error("Result for workunit %s was not"
                         " processed by any task."
                         " Setting it to status CANCELLED",
                         message.get_wu_id())
            self.wuar.cancel(message.get_wu_id())
        return was_received


# One entry in the WU DB, including the text with the WU contents
# (FILEs, COMMANDs, etc.) and info about the progress on this WU (when and
# to whom assigned, received, etc.)

# - wuid is the unique wuid of the workunit
# - status is a status code as defined in WuStatus
# - data is the str containing the text of the workunit
# - timecreated is the string containing the date and time of when the
#   WU was added to the db
# - timeassigned is the ... of when the WU was assigned to a client
# - assignedclient is the clientid of the client to which the WU was assigned
# - timeresult is the ... of when a result for this WU was received
# - resultclient is the clientid of the client that uploaded a result
#   for this WU
# - errorcode is the exit status code of the first failed command, or 0
#   if none failed
# - timeverified is the ... of when the result was marked as verified


if __name__ == '__main__':
    import argparse

    queries = {"avail": ("Available workunits",
                         {"eq": {"status": WuStatus.AVAILABLE}}),
               "assigned": ("Assigned workunits",
                            {"eq": {"status": WuStatus.ASSIGNED}}),
               "receivedok": ("Received ok workunits",
                              {"eq": {"status": WuStatus.RECEIVED_OK}}),
               "receivederr": ("Received with error workunits",
                               {"eq": {"status": WuStatus.RECEIVED_ERROR}}),
               "verifiedok": ("Verified ok workunits",
                              {"eq": {"status": WuStatus.VERIFIED_OK}}),
               "verifiederr": ("Verified with error workunits",
                               {"eq": {"status": WuStatus.VERIFIED_ERROR}}),
               "cancelled": ("Cancelled workunits",
                             {"eq": {"status": WuStatus.CANCELLED}}),
               "all": ("All existing workunits", {})
               }

    parser = argparse.ArgumentParser()
    parser.add_argument('-dbfile', help='Name of the database file')
    parser.add_argument('-create', action="store_true",
                        help='Create the database tables if they do not exist')
    parser.add_argument('-add', action="store_true",
                        help='Add new workunits. Contents of WU(s) are '
                        'read from stdin, separated by blank line')
    parser.add_argument('-assign', nargs=1, metavar='clientid',
                        help='Assign an available WU to clientid')
    parser.add_argument('-cancel', action="store_true",
                        help='Cancel selected WUs')
    parser.add_argument('-expire', action="store_true",
                        help='Expire selected WUs')
    # parser.add_argument('-setstatus', metavar='STATUS',
    #     help='Forcibly set selected workunits to status (integer)')
    parser.add_argument('-prio', metavar='N',
                        help='If used with -add, newly added WUs '
                        'receive priority N')
    parser.add_argument('-limit', metavar='N',
                        help='Limit number of records in queries',
                        default=None)
    parser.add_argument('-result', nargs=6,
                        metavar=('wuid', 'clientid', 'filename', 'filepath',
                                 'filetype', 'command'),
                        help='Return a result for wu from client')
    parser.add_argument('-test', action="store_true",
                        help='Run some self tests')
    parser.add_argument('-debug', help='Set debugging level')
    parser.add_argument('-setdict', nargs=4,
                        metavar=("dictname", "keyname", "type", "keyvalue"),
                        help='Set an entry of a DB-backed dictionary')

    parser.add_argument('-wuid', help="Select workunit with given id",
                        metavar="WUID")
    for arg in queries:
        parser.add_argument('-' + arg, action="store_true", required=False,
                            help="Select %s" % queries[arg][0].lower())
    parser.add_argument('-dump', nargs='?', default=None, const="all",
                        metavar="FIELD",
                        help='Dump WU contents, optionally a single field')
    parser.add_argument('-sort', metavar="FIELD",
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
            print(count)
        else:
            wus = db_pool.query(limit=args["limit"], **condition)
            if wus is None:
                print("0")
            else:
                print(len(wus))
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
            db_pool.set_status(WuStatus.NEED_RESUBMIT,
                               commit=True, **condition)
        # if args["setstatus"]:
        #    db_pool.set_status(int(args["setstatus"]), **condition)

    # Dict manipulation
    if args["setdict"]:
        (name, keyname, itemtype, keyvalue) = args["setdict"]
        # Type-cast value to the specified type
        value = getattr(__builtins__, itemtype)(keyvalue)
        dbdict = DictDbAccess(dbname, name)
        dbdict[keyname] = value
        del dbdict

    # Functions for testing
    if args["assign"]:
        clientid = args["assign"][0]
        wus = db_pool.assign(clientid)

    if args["result"]:
        result = args["result"]
        db_pool.result(result.wuid, result.clientid, result[2:])

# Local Variables:
# version-control: t
# End:
