from cadofactor.database.base import READONLY, EXCLUSIVE
from cadofactor.database.base import conn_close
from cadofactor.database.table import DbTable
from cadofactor.database.factory import DBFactory

from collections.abc import MutableMapping
from collections.abc import Mapping


# The reason for the 300-character limit is that the sqrt_factors table
# contains the input number to be factored. As such, we must make sure
# that it's permitted to go at least as far as we intend to go. 200
# digits is definitely too small.
class DictDbTable(DbTable):
    fields = (("rowid", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"),
              ("kkey", "VARCHAR(300)", "UNIQUE NOT NULL"),
              ("type", "INTEGER", "NOT NULL"),
              ("value", "TEXT", "")
              )
    primarykey = fields[0][0]
    references = None

    def __init__(self, *args, name=None, **kwargs):
        self.tablename = name
        # index creation now always prepends the table name,
        # and appends "index"
        self.index = {"dictdb_kkey": ("kkey",)}  # useful ?
        super().__init__(*args, **kwargs)


class DictDbDirectAccess(MutableMapping):
    """ A DB-backed flat dictionary.

    Flat means that the value of each dictionary entry must be a type that
    the underlying DB understands, like integers, strings, etc., but not
    collections or other complex types.

    >>> conn = DBFactory('db:sqlite3://:memory:').connect()
    >>> d = DictDbDirectAccess(conn, 'test')
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
    >>> del d
    >>> d = DictDbDirectAccess(conn, 'test')
    >>> d == {'a': 2, 'b': '3'}
    True
    >>> del d['b']
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
    >>> del d
    >>> d = DictDbAccess(conn, 'test')
    >>> d == {'a': '3', 'b': 3.0, 'c': '4', 'd': True}
    True
    >>> d.clear(['a', 'd'])
    >>> d == {'b': 3.0, 'c': '4'}
    True
    >>> del d
    >>> d = DictDbAccess(conn, 'test')
    >>> d == {'b': 3.0, 'c': '4'}
    True
    >>> d.clear()
    >>> d == {}
    True
    >>> del d
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
        self._table = DictDbTable(name=name)
        # Create an empty table if none exists
        self._conn.harness_transaction(EXCLUSIVE, self._table.create)

    def get_cursor(self):
        return self._conn.cursor()

    def __getitem__(self, key):
        r, = self._conn.harness_transaction(READONLY,
                                            self._table.where,
                                            # cursor is implicitly added
                                            limit=1,
                                            eq=dict(kkey=key))
        return self.__convert_value(r)

    def _iter_raw(self):
        n = len(self)
        batch_size = 128
        for i in range(1, n + 1, batch_size):
            rows = self._conn.harness_transaction(READONLY,
                                                  self._table.where,
                                                  # cursor is implicitly added
                                                  limit=batch_size,
                                                  offset=i)
            for r in rows:
                yield (r["kkey"], self.__convert_value(r))

    def __iter__(self):
        for r, v in self._iter_raw():
            yield r

    def items(self):
        return self._iter_raw()

    def __len__(self):
        return self._conn.harness_transaction(READONLY, self._table.count)

    def __contains__(self, key):
        return bool(self._conn.harness_transaction(READONLY,
                                                   self._table.count,
                                                   # cursor is implicitly added
                                                   eq=dict(kkey=key)))

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

    # This is guaranteed to hit the sql table and no other method
    def _backend_getall(self):
        """ Reads the whole table and returns it as a dict """
        rows = self._conn.harness_transaction(READONLY, self._table.where)
        return {r["kkey"]: self.__convert_value(r) for r in rows}

    def __setitem_nocommit(self, cursor, key, value):
        """ Set dictionary key to value and update/insert into table,
        but don't commit. Cursor must be given
        """
        # Insert a new row
        # if self._table.count(cursor, eq=dict(kkey=key)):
        r = self._table.where(cursor, limit=1, eq=dict(kkey=key))
        if r:
            self._table.update(cursor,
                               dict(value=str(value),
                                    type=self.__get_type_idx(value)),
                               eq=dict(kkey=key))
        else:
            self._table.insert(cursor,
                               dict(kkey=key,
                                    value=str(value),
                                    type=self.__get_type_idx(value)))

    def __setitem__(self, key, value):
        """ Access by indexing, e.g., d["foo"]. Always commits """
        self._conn.harness_transaction(EXCLUSIVE,
                                       self.__setitem_nocommit,
                                       # cursor is implicitly added
                                       key, value)

    def __delitem__(self, key, commit=True):
        """ Delete a key from the dictionary """
        self._conn.harness_transaction(EXCLUSIVE,
                                       self._table.delete,
                                       eq=dict(kkey=key))

    def setdefault(self, key, default=None, commit=True):
        ''' Setdefault function that allows a mapping as input

        Values from default dict are merged into self, *not* overwriting
        existing values in self '''
        if key is None and isinstance(default, Mapping):
            update = {key: default[key]
                      for key in default
                      if key not in self}
            if update:
                self.update(update, commit=commit)
            return None
        elif key not in self:
            self.update({key: default}, commit=commit)
        return self[key]

    def update(self, other, commit=True):
        def t(cursor):
            for (key, value) in other.items():
                self.__setitem_nocommit(cursor, key, value)
        self._conn.harness_transaction(EXCLUSIVE, t)

    def clear(self, args=None, commit=True):
        """ Overridden clear that allows removing several keys atomically """
        def t(cursor):
            if args is None:
                self._table.delete(cursor)
            else:
                for key in args:
                    self._table.delete(cursor, eq=dict(kkey=key))

        self._conn.harness_transaction(EXCLUSIVE, t)


class DictDbCachedAccess(DictDbDirectAccess):
    """
    This is the same as DictDbDirectAccess, except that everything is
    stored in a cached dictionary:
    a copy of all the data in the table is kept in memory; read accesses
    are always served from the in-memory dict. Write accesses write
    through to the DB.
    """

    types = (str, int, float, bool)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._data = self._backend_getall()

    # Implement the abstract methods defined by
    # collections.abc.MutableMapping All but __del__ and __setitem__ are
    # simply passed through to the self._data dictionary

    def __getitem__(self, key):
        return self._data.__getitem__(key)

    def __iter__(self):
        return self._data.__iter__()

    def __len__(self):
        return self._data.__len__()

    def __str__(self):
        return self._data.__str__()

    def __contains__(self, key):
        return key in self._data

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        # Update the in-memory dict
        self._data[key] = value

    def __delitem__(self, key, commit=True):
        super().__delitem__(key)
        # Update the in-memory dict
        del self._data[key]

    def setdefault(self, key, default=None, commit=True):
        '''
        note that key=None is a house addition
        '''
        if key is None and isinstance(default, Mapping):
            super().setdefault(key=key, default=default, commit=commit)
            return
        elif key not in self:
            self._data[key] = super().setdefault(key=key,
                                                 default=default,
                                                 commit=commit)
        return self._data[key]

    def update(self, other, commit=True):
        super().update(other, commit=commit)
        for (key, value) in other.items():
            self._data[key] = value

    def clear(self, args=None, commit=True):
        """ Overridden clear that allows removing several keys atomically """
        super().clear(args=args, commit=commit)
        if args is None:
            self._data.clear()
        else:
            for key in args:
                del self._data[key]


DictDbAccess = DictDbCachedAccess
