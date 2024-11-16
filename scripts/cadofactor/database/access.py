from cadofactor.database.dictdb import DictDbAccess


class IdMap(object):
    """
    Ensures that DB-backed dictionaries of the same table name are
    instantiated only once. This is made so that multiple control flows
    _within the same thread_ can have a one-stop shop to access this
    database table. It is _not_ wise to share this among threads, since
    it will crash!

    Since we assume that we have only one thread accessing the database,
    it is ok to use DictDbAccess (a.k.a. DictDbCachedAccess). But if we
    want to go multithread, which is conceivable, we would need
    DictDbDirectAccess, at some performance cost.
    """
    def __init__(self):
        self.db_dicts = {}

    def make_db_dict(self, db, name):
        key = name
        if key not in self.db_dicts:
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
