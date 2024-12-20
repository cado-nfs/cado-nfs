import sys
import sqlite3
from cadofactor.database.base import DB_base
from cadofactor.database.base import CursorWrapperBase
from cadofactor.database.base import ConnectionWrapperBase
from cadofactor.database.base import pending_transactions
from cadofactor.database.base import TransactionAborted
from cadofactor.database.base import logger


class DB_SQLite(DB_base):
    class CursorWrapper(CursorWrapperBase):
        @property
        def cursor(self):
            return self.__cursor

        @property
        def connection(self):
            return self.cursor.connection

        def __init__(self, cursor, *args, **kwargs):
            self.__cursor = cursor
            super().__init__(*args, **kwargs)

        def try_catch_execute(self, command, values):
            try:
                super().try_catch_execute(command, values)
            except (sqlite3.OperationalError, sqlite3.DatabaseError) as e:
                if str(e) == "database disk image is malformed" or \
                   str(e) == "disk I/O error":
                    logger.critical("sqlite3 reports errors"
                                    " accessing the database.")
                    logger.critical("Database file may have gotten corrupted,"
                                    " or maybe filesystem does not properly"
                                    " support  file locking.")
                    raise
                if str(e) != "database is locked":
                    raise
                raise TransactionAborted(command, values)

    class ConnectionWrapper(sqlite3.Connection, ConnectionWrapperBase):
        def cursor(self):
            return DB_SQLite.CursorWrapper(super().cursor())

        # https://docs.python.org/3/library/sqlite3.html#transaction-control
        def __init__(self, path, *args, **kwargs):
            ac = {}
            if sys.hexversion >= 0x030c0000:
                ac['autocommit'] = sqlite3.LEGACY_TRANSACTION_CONTROL
            super().__init__(path, *args,
                             # isolation_level=None,
                             **ac,
                             **kwargs)
            self.pending = pending_transactions.new_db(self)

    def connect(self):
        c = self.ConnectionWrapper(self.path)
        self.advertise_connection()
        return c

    # FIXME I think that in any case the sqlite3 module ends up creating
    # the db, no ?
    def __init__(self, uri, create=False):
        super().__init__(uri, backend_pattern="sqlite3?")
        self.path = self.db_name
