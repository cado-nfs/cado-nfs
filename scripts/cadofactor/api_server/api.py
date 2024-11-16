import sys
import flask
import logging
import os.path
import json
import time
import re
import socket
from ipaddress import ip_address, ip_network
import threading

from cadofactor import wudb
from werkzeug.utils import secure_filename
from werkzeug.exceptions import HTTPException

logger = logging.getLogger(__name__)

# some resources:
# https://medium.com/@DanKaranja/building-api-documentation-in-flask-with-swagger-a-step-by-step-guide-59a453509e2f

try:
    from flasgger import Swagger
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


class AppWrapper(flask.Flask):
    """
    This object embeds a flask.Flask app, and adds the minimal layer of
    things that we need for interaction with the rest.
    """

    def __init__(self, app_name):
        super().__init__(app_name)
        if 'flasgger' in sys.modules:
            self.swagger = Swagger(self, template=swagger_template)
        else:
            self.swagger = None
        self.database_uri = None
        self._db_connection_pool = {}
        self.serving_wus = None
        self.no_work_available = None
        self.only_registered = None
        self.upload_dir_provider = None
        self.upload_cycle = None
        self.timeout_hint = None
        self.wdir = None
        self.whitelist = None

    def hook_to_rest_of_code(self,
                             database_uri,
                             upload_dir_provider,
                             timeout_hint,
                             whitelist=None,
                             only_registered=True):
        self.database_uri = database_uri
        self._db_connection_pool = {}
        self.wuaccess = None
        self.serving_wus = True
        self.no_work_available = False
        self.only_registered = only_registered
        self.upload_dir_provider = upload_dir_provider
        self.upload_cycle = 0
        self.timeout_hint = timeout_hint
        self.whitelist = []
        for w in whitelist:
            try:
                self.whitelist.append(f"{socket.gethostbyname(w)}/32")
            except socket.gaierror:
                self.whitelist.append(w)
        logger.info(f"server whitelist is {self.whitelist}")
        self.wdir = None

    def _get_db_things(self):
        tid = threading.current_thread().ident
        m = self._db_connection_pool.get(tid)
        if m is None:
            c = self.database_uri.connect()
            wuaccess = wudb.WuAccess(c)
            files = wudb.DictDbDirectAccess(c, 'server_registered_filenames')
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


app = AppWrapper(__name__)


@app.errorhandler(HTTPException)
def handle_exception(e):
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
    return response


@app.before_request
def limit_remote_addr():
    """
    Implements ip filtering
    """
    peer = flask.request.remote_addr
    if app.whitelist is None:
        logger.error(f'blocked incoming request from {peer} ;'
                     ' NOTE: no whitelist is configured, all ip addresses'
                     ' are blocked anyway. You probably want to add'
                     ' a server.whitelist= argument')
        flask.abort(403)
    for net in app.whitelist:
        if ip_address(peer) in ip_network(net):
            return
    logger.error(f'blocked incoming request from {peer}')
    flask.abort(403)  # Forbidden


@app.route("/")
def hello_world():
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

    return {'message': "Hello, World!"}


@app.route("/db")
def get_db_name():
    """
    This returns misc stuff, mostly for debugging.
    ---
    tags:
      - misc
    produces:
      - application/json
    responses:
      200:
        description: A successful response
    """

    d = wudb.DictDbAccess(app.get_db_connection(), 'tasks')
    return {'message': app.database_uri, 'data': d._getall()}


@app.route("/workunit")
def getwu():
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
    if not (clientid := flask.request.form.get('clientid')):
        logger.debug(f"got {flask.request.path}"
                     f" from {flask.request.remote_addr}"
                     f" (no client identification provided)")
        flask.abort(403, 'clientid must be provided')

    logger.debug(f"got {flask.request.path}"
                 f" from {flask.request.remote_addr}"
                 f" with client identification {clientid}")

    if not app.serving_wus:
        flask.abort(410, "Distributed computation finished")

    # we might want to make the timeout dependent on the clientid.
    wu = app.get_wuaccess().assign(clientid,
                                   timeout_hint=app.timeout_hint)

    if not wu:
        # This flag is to downgrade the logging level. Ugly.
        app.no_work_available = True
        flask.abort(404, "No work available")

    logger.info(f"Sending workunit {wu.get_id()}"
                f" to client {clientid}")

    # change the timeout to a deadline for the json that we return to the
    # client. Inside the database, we really only maintain a timeout
    # value.
    if 'timeout' in wu:
        wu['deadline'] = time.time() + float(wu['timeout'])
        del wu['timeout']

    return wu


# the old server has a (presumably never used) mechanism to forward a
# request to the database. It's doing fairly ad hoc stuff that rewrites
# some on-the-spot DSL into an SQL query. Clearly low priority, I'm not
# even sure I see a use for this.


@app.route('/file/<path:path>')
def get_file(path):
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
    d = app.get_registered_filenames()
    if d is None:
        if re.match('^/', path):
            full_path = path
        else:
            full_path = os.path.join(app.get_workdir(), path)
        logger.info(f"got request for arbitrary file {path},"
                    f" which resolves to {full_path}")
        if os.path.isfile(full_path):
            return flask.send_file(full_path)
    elif (file := d.get(path)) is not None:
        logger.info(f'got request for file {path}, which resolves to {file}')
        print(file)
        if os.path.isfile(file):
            dirname, basename = os.path.split(file)
            return flask.send_from_directory(dirname, basename)
    flask.abort(404, 'File not found')


@app.route('/files')
def get_all_files():
    """
    List the server files that are available for download
    ---
    """
    d = app.get_registered_filenames()
    if d is None:
        return {}
    else:
        return d._getall()


@app.route('/failure')
def report_failure():
    # who = flask.request.form['clientid']
    # wuid = flask.request.form['WUid']
    # errorcode = flask.request.form['errorcode']
    # failedcommand = flask.request.form['failedcommand']
    pass


@app.route('/upload', methods=['POST'])
def upload_file():
    clientid = flask.request.form.get('clientid')
    wuid = flask.request.form.get('WUid')
    errorcode = flask.request.form.get('errorcode')
    failedcommand = flask.request.form.get('failedcommand')

    if clientid is None or wuid is None:
        flask.abort(400, "missing WUid and/or clientid")

    # XXX uploaded_files should be a list of tuples.
    # "filename", "path", "type", "command"

    fileinfo = json.loads(flask.request.form.get('fileinfo', "{}"))
    logger.info("file info: %s", fileinfo)

    uploaded_files = []
    for fkey, f in flask.request.files.items():
        filename = secure_filename(f.filename)
        path = os.path.join(app.get_upload_folder(), filename)
        if (fi := fileinfo.get(filename)) is None:
            flask.abort(400, f"missing fileinfo for file {filename}")

        if os.path.isfile(path):
            logger.error(f"denied request from {clientid}"
                         f" to overwrite existing file {path}")
            flask.abort(403, "File already exists")
        f.save(path)
        tup = [filename, path, fi["key"]]
        if (c := fi.get('command')) is not None:
            tup.append(c)
        uploaded_files.append(tuple(tup))
        # logger.info(f"saved result file {path}"
        #             f" for {wuid} (source={clientid})")

    try:
        app.get_wuaccess().result(wuid,
                                  clientid,
                                  uploaded_files,
                                  errorcode,
                                  failedcommand)
    except wudb.StatusUpdateError:
        logger.warning(f'Workunit {wuid} was not currently assigned')
    else:
        logger.info(f'Workunit {wuid} completed,'
                    f' uploaded files {uploaded_files}')

    return {'message': 'upload completed'}
