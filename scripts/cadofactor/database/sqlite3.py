import sys
import sqlite3
from cadofactor.database.base import DB_base
from cadofactor.database.base import CursorWrapperBase
from cadofactor.database.base import TransactionWrapper
from cadofactor.database.base import pending_transactions


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

    class ConnectionWrapper(sqlite3.Connection):
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

        def transaction(self, mode=None):
            return TransactionWrapper(self.cursor(), mode=mode)

    def connect(self):
        c = self.ConnectionWrapper(self.path)
        self.advertise_connection()
        return c

    # FIXME I think that in any case the sqlite3 module ends up creating
    # the db, no ?
    def __init__(self, uri, create=False):
        super().__init__(uri, backend_pattern="sqlite3?")
        self.path = self.db_name
