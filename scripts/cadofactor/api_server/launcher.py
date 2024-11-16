import logging
import multiprocessing
import sys
from cadofactor.cadofactor_tools import UploadDirProvider
from cadofactor import cadologger
from werkzeug import serving
from . import api
from .api import app

try:
    import gunicorn.app.base
except ImportError:
    pass

try:
    # interestingly, it seems that I don't even need the ssl module at all.
    import ssl  # noqa: F401
except ModuleNotFoundError:
    pass

HAVE_SSL = 'ssl' in sys.modules

if HAVE_SSL:
    from cadofactor.cadofactor_tools.certificate import \
            get_server_alternate_names, \
            create_certificate, \
            get_certificate_hash

# This implements the Cado-NFS api server. Some useful
# documentation sources that I perused in order to code this.
#
# https://blog.miguelgrinberg.com/post/running-your-flask-application-over-https


class ServerLauncher(object):
    """
    This is the cado-nfs api server. Its job is only to handle
    incoming requests, and communicate with the database. In turn, the
    rest of cado-nfs communicates with the database. (the wait loop is in
    cadotask.ClientServerTask.wait)
    """

    def __init__(self,
                 serveraddress,
                 serverport,
                 dbdata,
                 threaded=None,
                 debug=False,
                 uploaddir=None,
                 nrsubdir=None,
                 # scriptdir=None,
                 only_registered=True,
                 cafile=None,
                 whitelist=None,
                 timeout_hint=None,
                 # linger_before_quit=False
                 ):
        # some parameters are currently not targeted by the
        # implementation as it is now.

        """
        This starts the server at serveraddress:serverport, using the dbdata
        connection to the database.

        dbdata is a wudb.DBFactory object, which is basically a wrapper
        around an uri. Workers can use to make a fresh
        connection to the database by calling its connect() method. Note
        that dbdata itself is not a connection object!

        The threaded argument can be used to run a production WSGI server
        with gunicorn *if* this software is available. If not, we fall
        back to using a single-thread flask server, which should be
        sufficient for most purposes. Note that threaded=None is used as
        a default behaviour, where a threaded server is used
        opportunistically unless debug mode is on.

        The debug argument spawns a more noisy server.

        There are a few situations where it's desirable for the api server to
        process requests directly without even talking to the db. Namely, we
        want it to be able to return select files from the workdir (factor
        bases and such). Those are pulled from the database.
        The only_registered flag is a priori always equal to True. If False,
        a blanket access to all files is offered, which is probably unwise.

        TLS is used if cafile is given, in which case an ad hoc certificate
        is written to this file name.
        """

        self.name = "API server"
        self.logger = logging.getLogger(self.name)
        self.logger.setLevel(logging.NOTSET)

        self.address = serveraddress if serveraddress else "0.0.0.0"
        self.port = serverport

        self.cafile = cafile

        if not HAVE_SSL:
            if self.cafile:
                self.logger.warning("ssl not available,"
                                    f" cafile={cafile} ignored")
            self.cafile = None

        if self.cafile is not None:
            if 'ssl' not in sys.modules:
                raise RuntimeError("python ssl module is missing")
            SAN = get_server_alternate_names(serveraddress)
            self._ssl_context = create_certificate(self.cafile,
                                                   self.address,
                                                   SAN)
        else:
            self._ssl_context = None

        # https://stackoverflow.com/questions/70396641/how-to-run-gunicorn-inside-python-not-as-a-command-line
        if threaded is None:
            threaded = False if debug else True

        # inject fields to the api app object.

        app.hook_to_rest_of_code(dbdata,
                                 UploadDirProvider(uploaddir, nrsubdir),
                                 timeout_hint,
                                 whitelist=whitelist,
                                 only_registered=only_registered
                                 )

        api.logger.addHandler(cadologger.ScreenHandler(lvl=logging.NOTSET))

        if False and threaded and 'gunicorn' in sys.modules:
            class StandaloneApplication(gunicorn.app.base.BaseApplication):
                def __init__(self, app, options={}):
                    self.options = options
                    self.application = app
                    super().__init__()

                def load_config(self):
                    config = {key: value
                              for key, value in self.options.items()
                              if key in self.cfg.settings and value is not None
                              }
                    for key, value in config.items():
                        self.cfg.set(key.lower(), value)

                def load(self):
                    return self.application

            self.logger.info("Running from Gunicorn")

            if type(threaded) is int:
                nthreads = threaded
            else:
                nthreads = min(multiprocessing.cpu_count() * 2, 4) + 1

            def pre_fork(server, worker):
                print(f"pre-fork server {server} worker {worker}",
                      file=sys.stderr)

            options = {
                    'bind': f"{self.address}:{serverport}",
                    'workers': nthreads,
                    # 'threads': number_of_workers(),
                    'timeout': 120,
                    'pre_fork': pre_fork,
            }

            if self._ssl_context is not None:
                options['certfile'] = self._ssl_context[0]
                options['keyfile'] = self._ssl_context[1]

            self._run_object = StandaloneApplication(app, options)
            self._run_args = []
            self._run_kwargs = {}
        elif False:
            if threaded:
                self.logger.warning("Note: cannot run threaded server since"
                                    " we do not have gunicorn.")
            self.logger.info("Running from flask")
            self._run_object = app
            self._run_args = []
            self._run_kwargs = dict(host=serveraddress,
                                    port=serverport,
                                    debug=debug,
                                    use_reloader=False)
            if self._ssl_context is not None:
                self._run_kwargs['ssl_context'] = self
        else:
            class Foo(object):
                def __init__(self, address, port, app, **kwargs):
                    self.options = options
                    self.application = app
                    self.server = serving.make_server(host=address,
                                                      port=port,
                                                      app=self.application,
                                                      **kwargs)
                    self.port = self.server.socket.getsockname()[1]

                def load_config(self):
                    pass

                def load(self):
                    return self.application

                def run(self):
                    from threading import Thread
                    thread = Thread(target=self.server.serve_forever)
                    thread.daemon = True
                    thread.start()

            options = {}

            if threaded:
                if type(threaded) is int:
                    nthreads = threaded
                else:
                    nthreads = min(multiprocessing.cpu_count() * 2, 4) + 1

                options['threaded'] = nthreads

            if self._ssl_context is not None:
                options['ssl_context'] = self._ssl_context

            self._run_object = Foo(self.address, serverport, app, **options)
            self._run_args = []
            self._run_kwargs = {}
            self.logger.info("Running from werkzeug")

            scheme = 'https' if self._ssl_context else 'https'
            bound_address = self._run_object.server.socket.getsockname()[0]
            bound_port = self._run_object.server.socket.getsockname()[1]
            self.url = "%s://%s:%d" % (scheme, bound_address, bound_port)

            self.logger.info(f"Running on {self.url} (Press CTRL+C to quit))")
            self.port = bound_port
            self.url = self.url.replace('0.0.0.0', 'localhost')
            self.server = self._run_object.server

    def serve(self):
        self._run_object.run(*self._run_args, **self._run_kwargs)
        # document the connection procedure.
        connection_parameters = "--server=%s" % self.url
        if self._ssl_context:
            h = get_certificate_hash(self._ssl_context[0])
            connection_parameters += " --certsha1=%s" % h

        self.logger.info("You can start additional cado-nfs-client.py scripts"
                         " with parameters: " + connection_parameters)
        self.logger.info("If you want to start additional clients, remember "
                         "to add their hosts to server.whitelist")

    def get_cert_sha1(self):
        if self._ssl_context:
            return get_certificate_hash(self._ssl_context[0])
        else:
            return None

    def get_url(self, origin=None):
        # if origin=localhost, then we want to insist on using localhost
        # as the peer name. In fairness, I'm not sure it's really
        # important. It's only used in the all-clients-on-localhost case,
        # and then the default listen address of 0.0.0.0 gets rewritten
        # to localhost anyway.
        return self.url

    def get_port(self):
        return self.port

    def stop_serving_wus(self):
        self.logger.info("Got notification to stop serving Workunits")
        app.serving_wus = False

    def shutdown(self, exc=None):
        if exc is not None:
            self.logger.info("Shutting down server on exception %s", exc)
        else:
            self.logger.info("Shutting down server")
        # app.shutdown()
        self.server.shutdown()
