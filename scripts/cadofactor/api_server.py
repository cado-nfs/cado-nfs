"""
This implements the Cado-NFS api server. Some useful
documentation sources that I perused in order to code this.

https://blog.miguelgrinberg.com/post/running-your-flask-application-over-https
https://medium.com/@DanKaranja/building-api-documentation-in-flask-with-swagger-a-step-by-step-guide-59a453509e2f
https://stackoverflow.com/questions/39853643/is-it-possible-to-mark-a-method-with-a-route-in-flask

The new cado-nfs api server is built upon:

 - flask (python3-flask) in order to actually build an api
 - optionally, gunicorn (python3-gunicorn) in order to make a
   multithreaded server from it. (still WIP)
 - optionally, flasgger (python3-flasgger) in order to provide online api
   docs.
"""

import flask
import json
import logging
import multiprocessing
import os.path
import re
import socket
import sys
import threading
import time

from cadofactor.cadofactor_tools import UploadDirProvider
from cadofactor.database import DictDbDirectAccess
from cadofactor import wudb
from werkzeug import serving

from ipaddress import ip_address, ip_network

from werkzeug.utils import secure_filename
from werkzeug.exceptions import HTTPException


try:
    from flasgger import Swagger
except ModuleNotFoundError:
    pass

try:
    import gunicorn.app.base
except ImportError:
    pass

try:
    # interestingly, it seems that I don't even need the ssl module at all.
    import ssl  # noqa: F401
    from cadofactor.cadofactor_tools.certificate import \
        get_server_alternate_names, \
        create_certificate, \
        get_certificate_hash
except ModuleNotFoundError:
    pass

swagger_template = {
    "swagger": "2.0",
    "info": {
        "title": "Cado-NFS api",
        "description": "API for Cado-NFS server",
        "contact": {
            "name": "The Cado-NFS development team",
            "url": "https://cado-nfs.inria.fr",
            },
        "version": "0.1",
        "license": {
            "name": "GNU LGPL",
            "url":
                "https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html",
            }
        },
    }

HAVE_SSL = 'ssl' in sys.modules


class ApiServer(flask.Flask):
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
        super().__init__(self.name)
        self.logger = logging.getLogger(self.name)

        if True:
            # check for #30119
            import importlib
            import importlib.metadata

            def Version(c):
                return tuple([int(x) for x in c.split('.')])

            vw = Version(importlib.metadata.version('werkzeug'))
            v313 = Version("3.1.3")
            if vw <= v313:
                message = r"""
The flask package that is used here is affected by a bug
in the werkzeug package: See:
 - https://github.com/pallets/werkzeug/issues/3065
 - https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues/30119
It is expected that a future version of werkzeug (more
recent than 3.1.3) fixes this.
                """
                for m in message.strip().split("\n"):
                    self.logger.warning(m)

        self.address = serveraddress if serveraddress else "0.0.0.0"
        self.port = serverport

        self.cafile = cafile

        if not HAVE_SSL:
            if self.cafile:
                self.logger.warning("ssl not available,"
                                    f" cafile={cafile} ignored")
            self.cafile = None
            self._ssl_context = None

        if self.cafile is not None:
            if 'ssl' not in sys.modules:
                raise RuntimeError("python ssl module is missing")
            SAN = get_server_alternate_names(serveraddress)
            self._ssl_context = create_certificate(self.cafile,
                                                   self.address,
                                                   SAN)
            if self._ssl_context is None:
                self.logger.warning("ssl not available,"
                                    f" cafile={cafile} ignored")
                self.cafile = None
        else:
            self._ssl_context = None

        # https://stackoverflow.com/questions/70396641/how-to-run-gunicorn-inside-python-not-as-a-command-line
        if threaded is None:
            threaded = False if debug else True

        # inject fields to the api app object.

        if 'flasgger' in sys.modules:
            self.swagger = Swagger(self, template=swagger_template)
        else:
            self.swagger = None

        self.database_uri = dbdata
        self._db_connection_pool = {}
        self.wuaccess = None
        self.serving_wus = True
        self.no_work_available = False
        self.only_registered = only_registered
        self.upload_dir_provider = UploadDirProvider(uploaddir, nrsubdir)
        self.upload_cycle = 0
        self.timeout_hint = timeout_hint
        self.whitelist = []
        for w in whitelist or []:
            try:
                self.whitelist.append(f"{socket.gethostbyname(w)}/32")
            except socket.gaierror:
                self.whitelist.append(w)
        self.logger.info(f"server whitelist is {self.whitelist}")
        self.wdir = None

        self._route_endpoints()

        if False and threaded and 'gunicorn' in sys.modules:
            self._prepare_run_gunicorn(serveraddress, serverport, threaded)
        else:
            # self._prepare_run_flask_debug_server(
            #   serveraddress, serverport, threaded, debug=debug)
            self._prepare_run_werkzeug(serveraddress, serverport, threaded)

    def _prepare_run_gunicorn(self, serveraddress, serverport, threaded):
        class StandaloneApplication(gunicorn.app.base.BaseApplication):
            def __init__(self, app, options={}):
                self.options = options
                self.application = app
                super().__init__()
                # I'd love to be able to run with gunicorn, but
                # unfortunately I can't find a way to find the port it
                # binds to! The following is from the werkzeug server.
                # self.port = self.server.socket.getsockname()[1]

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
            # for debugging
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

        self._run_object = StandaloneApplication(self, options)
        self._run_args = []
        self._run_kwargs = {}

    def _prepare_run_flask_debug_server(self,
                                        serveraddress, serverport,
                                        threaded,
                                        debug=False):
        if threaded:
            self.logger.warning("Note: cannot run threaded server since"
                                " we do not have gunicorn.")
        self.logger.info("Running from flask")
        self._run_object = self
        self._run_args = []
        self._run_kwargs = dict(host=serveraddress,
                                port=serverport,
                                debug=debug,
                                use_reloader=False)
        if self._ssl_context is not None:
            self._run_kwargs['ssl_context'] = self

    def _prepare_run_werkzeug(self,
                              serveraddress, serverport,
                              threaded):
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

        nthreads = 1

        if threaded:
            if type(threaded) is int:
                nthreads = threaded
            else:
                nthreads = min(multiprocessing.cpu_count() * 2, 4) + 1

            options['threaded'] = nthreads

        if self._ssl_context is not None:
            options['ssl_context'] = self._ssl_context

        self._run_object = Foo(self.address, serverport, self, **options)
        self._run_args = []
        self._run_kwargs = {}
        self.logger.info("Running from werkzeug (%d thread(s))", nthreads)

        scheme = 'https' if self._ssl_context else 'http'
        bound_address = self._run_object.server.socket.getsockname()[0]
        bound_port = self._run_object.server.socket.getsockname()[1]
        self.url = "%s://%s:%d" % (scheme, bound_address, bound_port)

        self.port = bound_port
        self.url = self.url.replace('0.0.0.0', 'localhost')
        self.server = self._run_object.server

    def serve(self):
        self.logger.info(f"Running on {self.url} (Press CTRL+C to quit))")
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
        self.serving_wus = False

    def shutdown(self, exc=None):
        if exc is not None:
            self.logger.info("Shutting down server on exception %s", exc)
        else:
            self.logger.info("Shutting down server")
        # app.shutdown()
        self.server.shutdown()

    def _get_db_things(self):
        tid = threading.current_thread().ident
        m = self._db_connection_pool.get(tid)
        if m is None:
            c = self.database_uri.connect()
            wuaccess = wudb.WuAccess(c)
            files = DictDbDirectAccess(c, 'server_registered_filenames')
            self._db_connection_pool[tid] = (c, wuaccess, files)

        return self._db_connection_pool[tid]

    def get_db_connection(self):
        return self._get_db_things()[0]

    def get_wuaccess(self):
        return self._get_db_things()[1]

    def get_registered_filenames(self):
        return self._get_db_things()[2]

    def get_upload_folder(self):
        f = self.upload_dir_provider(self.upload_cycle)
        self.upload_cycle += 1
        return f

    def get_workdir(self):
        """
        Gets the main work directory as it is stored in the database
        """
        d = wudb.DictDbAccess(self.get_db_connection(), 'tasks')
        return d['workdir']

    def _route_endpoints(self):
        self.errorhandler(HTTPException)(self.api_errorhandler)
        self.before_request(self.api_limit_remote_addr)
        self.route("/")(self.api_hello_world)
        self.route("/workunit")(self.api_get_workunit)
        self.route("/file/<path:path>")(self.api_download_file)
        self.route("/files")(self.api_list_all_files)
        self.route("/upload", methods=["POST"])(self.api_upload_file)

    def api_errorhandler(self, e):
        """
        Return JSON instead of HTML for HTTP errors.
        """
        # start with the correct headers and status code from the error
        response = e.get_response()
        # replace the body with JSON
        response.data = json.dumps({
            "code": e.code,
            "name": e.name,
            "description": e.description,
        })
        response.content_type = "application/json"
        return response, 404

    def api_limit_remote_addr(self):
        """
        Implements ip filtering
        """
        peer = flask.request.remote_addr
        for net in self.whitelist or []:
            if ip_address(peer) in ip_network(net):
                return
        self.logger.error(f'blocked incoming request from {peer}')
        if self.whitelist is None:
            self.logger.error(' NOTE: no whitelist is configured,'
                              ' all ip addresses are blocked anyway.'
                              ' You probably want to add'
                              ' a server.whitelist= argument')
        flask.abort(403)  # Forbidden

    def api_hello_world(self):
        """
        This is an example endpoint that returns 'Hello, World!'
        ---
        tags:
          - misc
        responses:
            200:
                description: A successful response
                examples:
                    application/json: "Hello, World!"
        """

        resp = {'message': "Hello, World!"}

        return flask.json.jsonify(resp), 200

    def api_get_workunit(self):
        """
        This returns a fresh workunit.
        ---
        description: This returns a fresh workunit. The client must provide its
                     name with the clientid parameter. A 404 is return is the
                     server has not yet provisioned fresh workunits.
        parameters:
          - name: clientid
            in: query
            description: client-defined identifier
            required: true
            type: string
        produces:
          - application/json
        responses:
          200:
            description: A fresh workunit that the client can work on.
          404:
            description: No work available. This is a priori temporary, since
                         we are only waiting for the server to fill the
                         database with further workunits that we can hand
                         over to clients.
          410:
            description: The distributed computation phase is over. Clients
                         should terminate now.
        """
        clientid = flask.request.form.get('clientid')
        if clientid is None:
            # self.logger.debug(f"got {flask.request.path}"
            #                   f" from {flask.request.remote_addr}"
            #                   f" (no client identification provided)")
            flask.abort(403, 'clientid must be provided')

        # self.logger.debug(f"got {flask.request.path}"
        #                   f" from {flask.request.remote_addr}"
        #                   f" with client identification {clientid}")

        if not self.serving_wus:
            flask.abort(410, "Distributed computation finished")

        # we might want to make the timeout dependent on the clientid.
        wu = self.get_wuaccess().assign(clientid,
                                        timeout_hint=self.timeout_hint)

        if not wu:
            # This flag is to downgrade the logging level. Ugly.
            self.no_work_available = True
            flask.abort(404, "No work available")

        self.logger.info(f"Sending workunit {wu.get_id()}"
                         f" to client {clientid}")

        # change the timeout to a deadline for the json that we return to the
        # client. Inside the database, we really only maintain a timeout
        # value.
        if 'timeout' in wu:
            wu['deadline'] = time.time() + float(wu['timeout'])
            del wu['timeout']

        return flask.json.jsonify(wu), 200

    def api_download_file(self, path):
        """
        Pulls a file from the server.
        ---
        description: This pulls a file from the server, if this file is
                     registered for download
        produces:
          - application/octet-stream
        parameters:
          - name: path
            in: path
            description: file path as recognized by server
            required: true
            type: string
        """
        d = self.get_registered_filenames()
        if d is None:
            if re.match('^/', path):
                full_path = path
            else:
                full_path = os.path.join(self.get_workdir(), path)
            self.logger.info(f"got request for arbitrary file {path},"
                             f" which resolves to {full_path}")
            if os.path.isfile(full_path):
                return flask.send_file(full_path)
        else:
            file = d.get(path)
            if file is not None:
                self.logger.info(f"got request for file {path},"
                                 f" which resolves to {file}")
                if os.path.isfile(file):
                    dirname, basename = os.path.split(file)
                    return flask.send_from_directory(dirname, basename)
        flask.abort(404, 'File not found')

    def api_list_all_files(self):
        """
        List the server files that are available for download
        ---
        """
        d = self.get_registered_filenames()

        d = {} if d is None else dict(d)

        return flask.json.jsonify(d), 200

    def api_upload_file(self):
        clientid = flask.request.form.get('clientid')
        wuid = flask.request.form.get('WUid')
        errorcode = flask.request.form.get('errorcode')
        failedcommand = flask.request.form.get('failedcommand')

        if clientid is None or wuid is None:
            flask.abort(400, "missing WUid and/or clientid")

        # XXX uploaded_files should be a list of tuples.
        # "filename", "path", "type", "command"

        fileinfo = json.loads(flask.request.form.get('fileinfo', "{}"))

        uploaded_files = []
        for fkey, f in flask.request.files.items():
            filename = secure_filename(f.filename)
            path = os.path.join(self.get_upload_folder(), filename)
            fi = fileinfo.get(filename)
            if fi is None:
                self.logger.error("Incomplete answer from client:"
                                  " files=%s, data=%s",
                                  flask.request.files.keys(),
                                  json.dumps(fileinfo, indent=4))
                flask.abort(400, f"missing fileinfo for file {filename}")

            if os.path.isfile(path):
                self.logger.error(f"denied request from {clientid}"
                                  f" to overwrite existing file {path}")
                flask.abort(403, "File already exists")
            f.save(path)
            tup = [filename, path, fi["key"]]
            c = fi.get('command')
            if c is not None:
                tup.append(c)
            uploaded_files.append(tuple(tup))
            # self.logger.info(f"saved result file {path}"
            #                  f" for {wuid} (source={clientid})")

        try:
            self.get_wuaccess().result(wuid,
                                       clientid,
                                       uploaded_files,
                                       errorcode,
                                       failedcommand)
        except wudb.StatusUpdateError:
            self.logger.warning(f'Workunit {wuid} was not currently assigned')
        else:
            self.logger.debug(f'Workunit {wuid} completed,'
                              f' uploaded files {uploaded_files}')

        resp = {'message': 'upload completed'}

        return flask.json.jsonify(resp), 200
