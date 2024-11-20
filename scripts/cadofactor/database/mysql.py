import re
import copy
import mysql.connector

from cadofactor.database.base import DB_base
from cadofactor.database.base import CursorWrapperBase
from cadofactor.database.base import ConnectionWrapperBase
from cadofactor.database.base import pending_transactions

from cadofactor.database.base import logger
from cadofactor.database.base import TransactionAborted


class DB_MySQL(DB_base):
    class CursorWrapper(CursorWrapperBase):
        @property
        def parameter_auto_increment(self):
            return "%s"

        @property
        def _string_translations(self):
            return [
                    ('\\bASC\\b', "AUTO_INCREMENT"),
                    # create index if not exists seems to be okay with
                    # mariadb 11.x at least
                    # ('\\bCREATE INDEX IF NOT EXISTS\\b', "CREATE INDEX"),
                    ('\\bBEGIN EXCLUSIVE\\b', "START TRANSACTION"),
                    ('\\bpurge\\b', "purgetable"),
            ]

        @property
        def cursor(self):
            return self.__cursor

        @property
        def connection(self):
            return self._connection

        def upgrade_command(self, cmd):
            if type(cmd) is tuple:
                return tuple([self.upgrade_command(u) for u in cmd])
            elif type(cmd) is list:
                return [self.upgrade_command(u) for u in cmd]
            else:
                if re.search(r"^SELECT", cmd):
                    logger.transaction("Upgrading SELECT with FOR UPDATE")
                    cmd = re.sub(r';$', ' FOR UPDATE;', cmd)
                return cmd

        def __init__(self, cursor, connection=None, *args, **kwargs):
            self._connection = connection
            self.__cursor = cursor
            super().__init__(*args, **kwargs)

        def try_catch_execute(self, command, values):
            try:
                super().try_catch_execute(command, values)
            except mysql.connector.errors.InternalError as e:
                # we only want to raise our custom exceptions in cases
                # that we expect. Yes, the sql state is a __string__
                if e.sqlstate == '40001':
                    raise TransactionAborted(command, values)
                raise

    class ConnectionWrapper(ConnectionWrapperBase):
        def _reconnect_anonymous(self):
            self._conn = mysql.connector.connect(**self.db_connect_args)

        def _reconnect(self):
            self._conn = mysql.connector.connect(database=self.db_name,
                                                 **self.db_connect_args)
            cursor = self._conn.cursor()
            cursor.execute('SET AUTOCOMMIT=0;')
            self._conn.commit()
            cursor.execute('SET TRANSACTION ISOLATION LEVEL SERIALIZABLE;')
            self._conn.commit()

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
            self.db_name = self._db_factory.db_name
            self.db_connect_args = copy.copy(self._db_factory.db_connect_args)

            if self.db_connect_args.get('host') == '_unix_socket':
                m = re.match(r"^(.*)/([^/]*)$", self.db_name)
                if not m:
                    raise RuntimeError(
                        "mysql _unix_socket requires the unix socket path"
                        " to be passed as a prefix to the database name")
                del self.db_connect_args['host']
                self.db_connect_args['unix_socket'] = m.group(1)
                self.db_connect_args['password'] = None
                assert self.db_connect_args['user'] is not None
                self.db_name = m.group(2)

            if create:
                try:
                    self._reconnect()
                except mysql.connector.errors.ProgrammingError:
                    # need to create the database first. Do it by
                    # hand, with a connection which starts without a
                    # database name.
                    logger.info("Creating database %s" % self.db_name)
                    self._reconnect_anonymous()
                    cursor = self._conn.cursor()
                    cursor.execute("CREATE DATABASE %s;" % self.db_name)
                    cursor.execute("USE %s;" % self.db_name)
                    # cursor.execute("SET autocommit = 1")
                    self._conn.commit()
            else:
                self._reconnect()

            # several connections may have transactions running
            # concurrently, I believe. The backend will serialize them
            # self.pending = pending_transactions.new_db(db_factory.uri)
            self.pending = pending_transactions.new_db(self)

            logger.info(f"database: {self.db_name},"
                        f" connect args: {self.db_connect_args}")

        def rollback(self):
            self._conn.rollback()

        def close(self):
            self._conn.close()

        def commit(self):
            self._conn.commit()

        def in_transaction(self):
            return self._conn.in_transaction

    def connect(self, *args, **kwargs):
        return self.ConnectionWrapper(self, *args, **kwargs)

    def __init__(self, uri, create=False):
        super().__init__(uri, backend_pattern="mysql")
        self.path = None
        if create:
            conn = self.connect(create=True)
            conn.close()
