from cadofactor.database.base import logger


class DbTable(object):
    """ A class template defining access methods to a database table """

    @staticmethod
    def _subdict(d, ell):
        """
        Returns a dictionary of those key:value pairs of d for which
        the keys from ell are defined.
        """
        if d is None:
            return None
        return {k: d[k] for k in d.keys() if k in ell}

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
            r = self.references  # referenced table
            fk = (r.getpk(), "INTEGER", "REFERENCES %s ( %s ) "
                  % (r.getname(), r.getpk()))
            fields.append(fk)
        cursor.create_table(self.tablename, fields)
        if self.references:
            # We always create an index on the foreign key
            cursor.create_index(self.tablename + "_pkindex", self.tablename,
                                (fk[0], ))
        for indexname in self.index:
            # cursor.create_index(self.tablename + "_" + indexname,
            #                     self.tablename, self.index[indexname])
            try:
                cursor.create_index(f"{self.tablename}_{indexname}_index",
                                    self.tablename, self.index[indexname])
            except Exception as e:
                logger.warning(e)
                pass

    def insert(self, cursor, values, foreign=None):
        """
        Insert a new row into this table. The column:value pairs are
        specified key:value pairs of the dictionary d.  The database's
        row id for the new entry is stored in d[primarykey]
        """
        d = self.dictextract(values)
        assert self.primarykey not in d or d[self.primarykey] is None
        # If a foreign key is specified in foreign, add it to the column
        # that is marked as being a foreign key
        if foreign:
            r = self.references.primarykey
            assert r not in d or d[r] is None
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

    def count(self, cursor, **conditions):
        return cursor.count(self.tablename, **conditions)
