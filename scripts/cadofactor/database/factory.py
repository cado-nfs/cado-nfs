from cadofactor.database.base import DB_base


class DBFactory(object):
    # This class initializes the database from the supplied db uri.
    # db:engine:[//[user[:password]@][host][:port]/][dbname][?params][#fragment]
    def __init__(self, uri, *args, **kwargs):
        self.uri = uri
        self.base = None
        error = {}
        sc = DB_base.__subclasses__()
        for c in sc:
            # logger.info("Trying database back-end %s (among %d)"
            #             % (c, len(sc)))
            try:
                self.base = c(uri, *args, **kwargs)
                break
            except ValueError as err:
                error[str(c)] = err
                pass
        if self.base is None:
            msg = "Cannot use database URI %s" % uri
            msg += "\n" + "Messages received from %d backends:" % len(sc)
            for c in error.keys():
                msg += "\n" + "Error from %s: %s" % (c, error[c])
            raise ValueError(msg)

    def connect(self):
        return self.base.connect()

    @property
    def uri_without_credentials(self):
        return self.base.uri_without_credentials

    @property
    def path(self):
        # TODO: remove
        return self.base.path
