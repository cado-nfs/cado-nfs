import mysql.connector

from cadofactor.database.base import DB_base
from cadofactor.database.base import CursorWrapperBase
from cadofactor.database.base import TransactionWrapper
from cadofactor.database.base import pending_transactions

from cadofactor.database.base import logger


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
            self.__cursor = cursor
            super().__init__(*args, **kwargs)

    class ConnectionWrapper(object):
        def _reconnect_anonymous(self):
            self._conn = mysql.connector.connect(
                **self._db_factory.db_connect_args)

        def _reconnect(self):
            self._conn = mysql.connector.connect(
                database=self._db_factory.db_name,
                **self._db_factory.db_connect_args)

        def cursor(self):
            # provide some retry capability. This must be done on the
            # connection object, since reconnecting changes the
            # connection member.
            for i in range(10):
                try:
                    c = self._conn.cursor()
                    break
                except mysql.connector.errors.OperationalError:
                    logger.warning("Got exception connecting"
                                   " to the database, retrying (#%d)" % i)
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
            self.pending = pending_transactions.new_db(db_factory.uri)

        def rollback(self):
            self._conn.rollback()

        def close(self):
            self._conn.close()

        def commit(self):
            self._conn.commit()

        def in_transaction(self):
            return self._conn.in_transaction

        def transaction(self, mode=None):
            return TransactionWrapper(self.cursor(), mode=mode)

    def connect(self, *args, **kwargs):
        return self.ConnectionWrapper(self, *args, **kwargs)

    def __init__(self, uri, create=False):
        super().__init__(uri, backend_pattern="mysql")
        self.path = None
        if create:
            conn = self.connect(create=True)
            conn.close()
