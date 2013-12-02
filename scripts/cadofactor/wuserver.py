#!/usr/bin/env python3

import http.server
import socket
import socketserver
from socket import error as socket_error
import os
import sys
import re
import io
import logging
import urllib.parse
import copy
import struct
from workunit import Workunit
import datetime
import wudb
import upload
import select
import errno
from subprocess import check_output, CalledProcessError, STDOUT
try:
    import ssl
    HAVE_SSL = True
except ImportError:
    HAVE_SSL = False
    

# Import upload to get the shell environment variable name in which we should
# store the path to the upload directory

class FixedHTTPServer(http.server.HTTPServer):
    """ Work-around class for http.server.HTTPServer that handles EINTR """

    def __init__(self, *args, whitelist=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger = logging.getLogger("HTTP server")
        if whitelist is None:
            self.whitelist = None
        else:
            self.logger.info("Using whitelist: %s", ",".join(whitelist))
            self.whitelist = []
            for iprange in whitelist:
                mask = self.ipmask(iprange)
                if not mask:
                    raise ValueError("%s it not a valid IP range (must be "
                            "CIDR notation)" % iprange)
                self.whitelist.append(mask)

    @staticmethod
    def ipmask(iprange):
        """ Convert CIDR network range strings to a network address and a
        bit mask

        Convert text strings in the usual network mask form, e.g.,
        "192.168.3.0/24" to (netaddr, andmask) pairs which we can use to
        match IP addresses, i.e., the IP address is in the range iff
        (ip & andmask) == netaddr

        >>> FixedHTTPServer.ipmask("192.168.3.0/24")
        (3232236288, 4294967040)
        >>> FixedHTTPServer.ipmask("192.168.3.1")
        (3232236289, 4294967295)
        >>> FixedHTTPServer.ipmask("1.0.0.0/0")
        (0, 0)
        >>> FixedHTTPServer.ipmask("1.0.0.1/32")
        (16777217, 4294967295)
        >>> FixedHTTPServer.ipmask("localhost")
        (2130706433, 4294967295)
        >>> FixedHTTPServer.ipmask("1.0.0.256")
        >>> FixedHTTPServer.ipmask("1.0.0.0/33")
        """
        addr = iprange.split('/')
        if len(addr) < 2:
            addr.append(32)
        # Maybe it is a hostname. Try to resolve it
        try:
            addr[0] = socket.gethostbyname(addr[0])
        except socket.gaierror:
            return None
        try:
            packedIP = socket.inet_aton(addr[0])
        except socket.error:
            return None
        try:
            addr[1] = int(addr[1])
        except ValueError:
            return None
        if not 0 <= addr[1] <= 32:
            return None
        netaddr = struct.unpack("!L", packedIP)[0]
        andmask = 2**32 - 2**(32-addr[1])
        return (netaddr & andmask, andmask)

    def serve_forever(self, *args, **kwargs):
        """ Wrapper around http.server.HTTPServer.serve_forever() that
        restarts in case of EINTR.

        See http://bugs.python.org/issue7978
        """
        while True:
            try:
                return super().serve_forever(*args, **kwargs)
            except OSError as e:
                if e.errno != errno.EINTR:
                    raise
            except select.error as e:
                if e.args[0] != errno.EINTR:
                    raise

    def verify_request(self, request, client_address):
        """ Tests if the client's IP address is whitelisted

        If no whitelist is defined, always denies.
        """
        if not self.whitelist is None:
            # Use ipmask() to convert dotted string form of address to integer
            addr = self.ipmask(client_address[0])[0]
            for iprange in self.whitelist:
                if addr & iprange[1] == iprange[0]:
                    return True
        self.logger.warning("Connection from IP address %s rejected - "
                "not in server.whitelist", client_address[0])
        return False


class ThreadedHTTPServer(socketserver.ThreadingMixIn, FixedHTTPServer):
    """Handle requests in a separate thread."""


if HAVE_SSL:
    class HTTPSServer(FixedHTTPServer):
        def __init__(self, server_address, HandlerClass, *args, certfile=None,
                    keyfile=None, **kwargs):
            # Let TCPServer.__init__() call BaseServer.__init__() and create 
            # a self.socket attribute, but not bind and activate the socket
            super().__init__(server_address, HandlerClass, *args, 
                    bind_and_activate=False, **kwargs)
            # Create an SSL wrapper around the network socket in self.socket
            # First init an SSL context with the key
            ctx = ssl.SSLContext(ssl.PROTOCOL_SSLv23)
            ctx.load_cert_chain(certfile=certfile, keyfile=keyfile)
            # Now replace our orginal network socket with the SSL wrapper
            self.socket = ctx.wrap_socket(self.socket, server_side=True)
            # Finally bind and activate the SSLWrapper socket
            self.server_bind()
            self.server_activate()

    class ThreadedHTTPSServer(socketserver.ThreadingMixIn, HTTPSServer):
        """Handle requests in a separate thread."""

    BUGGY_SSLSOCKET_VERSIONS = [
        (3,2,0), (3,2,1), (3,2,2), (3,2,3),
        (3,3,0)
    ]
    if sys.version_info[0:3] in BUGGY_SSLSOCKET_VERSIONS:
        class FixedSSLSocket(ssl.SSLSocket):
            """ Wrapper class that applies the patch for issue 16357
            
            See http://bugs.python.org/issue16357
            """
            def accept(self):
                newsock, addr = socket.socket.accept(self)
                return (ssl.SSLSocket(sock=newsock,
                                      server_side=True,
                                      do_handshake_on_connect=
                                          self.do_handshake_on_connect,
                                      _context=self.context),
                        addr)
        ssl.SSLSocket = FixedSSLSocket


class HtmlGen(io.BytesIO):
    def __init__(self, encoding = None):
        super().__init__()
        if encoding is None:
            self.encoding = 'utf-8'
        else:
            self.encoding = encoding

    def header(self):
        self.write(
            b'<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" ' + 
            b'"http://www.w3.org/TR/html4/strict.dtd">\n' + 
            b'<html>\n' + 
            b'<head>\n' + 
            b'<meta http-equiv="content-type" content="text/html; ' + 
              b'charset=' + self.encoding.encode("ascii") + b'">\n' 
            b'<title>List of workunits</title>\n' + 
            b'</head>\n' + 
            b'<body>')

    def finish(self):
        self.write(b'</body>')

    def __bytes__(self):
        return self.getvalue()

    def get_len(self):
        return len(self.getvalue())

    def append(self, str):
        self.write(str.encode(self.encoding))

    def start_table(self, fields):
        self.append('<table border="1">\n<tr>')
        for h in fields:
            self.append('<th>' + h + '</th>')
        self.append('</tr>\n')

    def add_table_row(self, row):
        self.append('<tr>')
        for d in row:
            self.append('<td>' + str(d) + '</td>')
        self.append('</tr>\n')

    def end_table(self):
        self.append('</table>\n')

    def wu_row(self, wu, fields, cwd):
        arr = []
        for k in fields:
            if k == "files" and not wu["files"] is None:
                s = ""
                for f in wu["files"]:
                    path = f["path"]
                    if path.startswith(cwd):
                        path = path[len(cwd):]
                    s = s + '<a href="' + path + '">' + f["filename"] + \
                    '</a><br>'
                arr.append(s)
            else:
                arr.append(wu[k])
        self.add_table_row(arr)


class MyHandler(http.server.CGIHTTPRequestHandler):
    clientid_pattern = re.compile("^clientid=([\w.-]*)$")

    # The instance variable rbufsize is used by StreamRequestHandler to set
    # the buffer size of the read file object attached to the socket.
    # CGIHTTPRequestHandler sets this to 0, i.e., unbuffered, as it has to
    # pass the underlying file descriptor to a subprocess, so any data left
    # in the buffer would not be passed on to the subprocess.
    # We don't use subprocesses, thus we restore the original default 
    # buffering mode to avoid a huge performance hit.
    rbufsize=-1
    
    # Check that urlsplit() does not collapse paths which could prevent
    # registered_filenames from filtering path traversal attacks
    # See http://bugs.python.org/issue19435
    assert urllib.parse.urlsplit("http://foo//a").path == "//a"
    
    def __init__(self, *args, **kwargs):
        self.no_work_available = False
        super().__init__(*args, **kwargs)
    
    def log(self, lvl, format, *args, **kwargs):
        """ Interface to the logger class. 
            We add the client address (as a string) to the log record so the 
            logger can print that """
        # e = kwargs.copy()
        # e["address_string"] = self.address_string()
        format = '%s ' + format
        self.logger.log(lvl, format, self.address_string(), *args, **kwargs)

    # These three methods overwrite the corresponding methods from 
    # http.server.BaseHTTPRequestHandler
    # They just call self.log() with a numerical logging level added
    def log_message(self, format, *args, **kwargs):
        self.log(logging.INFO, format, *args, **kwargs)

    def log_request(self, code='-', size='-'):
        self.log(logging.DEBUG, '"%s" %s %s', self.requestline, code, size)

    def log_error(self, format, *args, **kwargs):
        # Log errors with WARNING level, except messages about no work being
        # available, as those are frequent and kinda spammy
        level = logging.DEBUG if self.no_work_available else logging.WARNING
        self.log(level, format, *args, **kwargs)

    def send_body(self, body):
        self.wfile.write(body)
        self.wfile.flush()
    
    def translate_path(self, path):
        """ Translate path in request URL to local file system, taking into 
        account registered file names.
        Overrides SimpleHTTPRequestHandler.translate_path(); paths that are not
        in registered_filenames are delegated to super().translate_path()
        """
        # Path in url always starts with '/'
        # print("translate_path(%s)" % path)
        relpath = path.lstrip('/')
        if relpath in self.registered_filenames:
            self.log(logging.DEBUG, "Translated path %s to %s", path, 
                     self.registered_filenames[relpath])
            return self.registered_filenames[relpath]
        self.log(logging.DEBUG, "Not translating path %s ", path)
        if self.only_registered:
            return None
        else:
            return super().translate_path(path)
    
    def do_GET(self):
        """Generates a workunit if request is cgi-bin/getwu, generates a status
        page is requested, otherwise calls parent class' do_GET()"""
        if self.is_cgi():
            if self.is_getwu():
                self.send_WU()
            elif self.is_getstatus():
                self.send_status()
            else:
                self.send_error(404, "GET for CGI scripts allowed only "
                                "for workunit or status page request")
        elif self.only_registered and \
                not self.path.lstrip('/') in self.registered_filenames:
                self.send_error(404, "Access restricted to registered file "
                                "names, %s is not registered" % self.path)
        else:
            try:
                super().do_GET()
            except socket_error as e:
                self.log_error("Connection error: %s", str(e))
        sys.stdout.flush()

    def do_POST(self):
        if self.is_upload():
            self.send_response(200, "Script output follows")
            # Python 3.3 needs flush_headers()
            try:
                self.flush_headers()
            except AttributeError:
                pass
            try:
                self.do_upload()
            except socket.error as e:
                self.log_error("%s", e)
        else:
            self.send_error(501, "POST request allowed only for uploads")
        sys.stdout.flush()
        sys.stderr.flush()

    def is_upload(self):
        """Test whether request is a file upload."""
        
        splitpath = urllib.parse.urlsplit(self.path)
        if self.command == 'POST' and self.is_cgi() and \
                splitpath.path == self.upload_path:
            return True
        return False

    def do_upload(self):
        if False:
            # This executes upload.py as a CGI script, however, that does
            # not work with SSL because the SSL-wrapper around the socket
            # is missing proper file descriptors for dup(), etc.
            # This does not work if rbufsize != 0.
            super().do_POST()
        else:
            # This uses the imported do_upload() function directly, without
            # spawning a subprocess, so that do_upload() can use the file-
            # like objects provided by the SSL wrapper. I should be faster,
            # too, by not having to execute a new Python interpreter and
            # script for each upload.
            # The CGI parsing code used in upload.py requires some shell
            # environment variables to be set according to the CGI
            # specificaton, such as CONTENT_LENGTH.
            # This is really slow if rbufsize == 0.
            env = self.create_env("upload.py")
            upload.do_upload(self.dbfilename, self.uploaddir, 
                    inputfp=self.rfile, output=self.wfile, environ=env)

    def create_env(self, scriptname, source_env=os.environ, query=None):
        """ Create a set of shell environment variables according to the CGI
        specification.
        
        Copied from the Python 3.2 http/server.py library file.
        """
        env = copy.deepcopy(source_env)
        env['SERVER_SOFTWARE'] = self.version_string()
        env['SERVER_NAME'] = self.server.server_name
        env['GATEWAY_INTERFACE'] = 'CGI/1.1'
        env['SERVER_PROTOCOL'] = self.protocol_version
        env['SERVER_PORT'] = str(self.server.server_port)
        env['REQUEST_METHOD'] = self.command
        uqrest = urllib.parse.unquote(scriptname)
        env['PATH_INFO'] = uqrest
        env['PATH_TRANSLATED'] = self.translate_path(uqrest)
        env['SCRIPT_NAME'] = scriptname
        if query:
            env['QUERY_STRING'] = query
        host = self.address_string()
        if host != self.client_address[0]:
            env['REMOTE_HOST'] = host
        env['REMOTE_ADDR'] = self.client_address[0]
        authorization = self.headers.get("authorization")
        if authorization:
            authorization = authorization.split()
            if len(authorization) == 2:
                import base64, binascii
                env['AUTH_TYPE'] = authorization[0]
                if authorization[0].lower() == "basic":
                    try:
                        authorization = authorization[1].encode('ascii')
                        authorization = base64.decodebytes(authorization).\
                                        decode('ascii')
                    except (binascii.Error, UnicodeError):
                        pass
                    else:
                        authorization = authorization.split(':')
                        if len(authorization) == 2:
                            env['REMOTE_USER'] = authorization[0]
        # XXX REMOTE_IDENT
        if self.headers.get('content-type') is None:
            env['CONTENT_TYPE'] = self.headers.get_content_type()
        else:
            env['CONTENT_TYPE'] = self.headers['content-type']
        length = self.headers.get('content-length')
        if length:
            env['CONTENT_LENGTH'] = length
        referer = self.headers.get('referer')
        if referer:
            env['HTTP_REFERER'] = referer
        accept = []
        for line in self.headers.getallmatchingheaders('accept'):
            if line[:1] in "\t\n\r ":
                accept.append(line.strip())
            else:
                accept = accept + line[7:].split(',')
        env['HTTP_ACCEPT'] = ','.join(accept)
        ua = self.headers.get('user-agent')
        if ua:
            env['HTTP_USER_AGENT'] = ua
        co = filter(None, self.headers.get_all('cookie', []))
        cookie_str = ', '.join(co)
        if cookie_str:
            env['HTTP_COOKIE'] = cookie_str
        # XXX Other HTTP_* headers
        # Since we're setting the env in the parent, provide empty
        # values to override previously set values
        for k in ('QUERY_STRING', 'REMOTE_HOST', 'CONTENT_LENGTH',
                  'HTTP_USER_AGENT', 'HTTP_COOKIE', 'HTTP_REFERER'):
            env.setdefault(k, "")
        return env


    def is_getwu(self):
        """Test whether request is for a new WU."""
        if not self.command == 'GET':
            return False
        splitpath = urllib.parse.urlsplit(self.path)
        # print("is_getwu(): path = %s, splitpath.path = %s" %
        #       (self.path, splitpath.path))
        return splitpath.path in [directory + '/getwu'
                                  for directory in self.cgi_directories]

    def is_getstatus(self):
        """Test whether request is for a a status page."""
        if not self.command == 'GET':
            return False
        splitpath = urllib.parse.urlsplit(self.path)
        # print("is_getwu(): path = %s, splitpath.path = %s" %
        #       (self.path, splitpath.path))
        return splitpath.path in [directory + '/status'
                                  for directory in self.cgi_directories]

    def guess_type(self, path):
        type = super().guess_type(path)
        # Use text/plain for files in upload, unless the type was properly 
        # identified
        # FIXME: testing the CWD here is wrong. If we want to expose files in
        # the upload directory, we either need to register all of them, or
        # register the upload directory and add directory name translation.
        cwd = os.getcwd().rstrip(os.sep)
        if type == "application/octet-stream" and \
                 path.startswith(cwd + os.sep + 'upload' + os.sep):
            return "text/plain"
        return type

    def send_WU(self):
        splitpath = urllib.parse.urlsplit(self.path)
        # print(str(splitpath))
        clientid_match = self.clientid_pattern.match(splitpath.query)
        if not clientid_match:
            return self.send_error(400, "No or malformed client id specified")
        clientid = clientid_match.group(1)
        
        if self.db_pool:
            wu_text = self.db_pool.assign(clientid)
        else:
            wu_text = wudb.WuAccess(self.dbfilename).assign(clientid)
        if not wu_text:
            # This flag is to downgrade the logging level. Ugly.
            self.no_work_available = True
            try:
                return self.send_error(404, "No work available")
            except socket_error as e:
                self.log_error("Connection error: %s", str(e))
                return
        
        self.log_message("Sending workunit " + Workunit(wu_text).get_id() + 
                         " to client " + clientid)
        # wu_text = wu.get_wu()
        try:
            self.send_response(200)
            self.send_header("Content-Type", "text/plain")
            self.send_header("Content-Length", len(wu_text))
            self.send_header("Cache-Control", "no-cache")
            self.end_headers()
            # FIXME: is ASCII enough for workunits? Is there any shell syntax
            # that needs more, or should we allow non-ASCII workunit names?
            self.send_body(bytes(wu_text, "ascii"))
        except socket_error as e:
            self.log_error("Connection error: %s", str(e))

    def send_status(self):
        self.send_query()
    
    def send_query(self):
        splitpath = urllib.parse.urlsplit(self.path)
        conditions = {}
        logging.debug("send_query(): query: %s", splitpath.query)
        if splitpath.query:
            # Parse query part into SELECT conditions
            # Now look at individual key=value pairs
            for query in splitpath.query.split("&"):
                query = urllib.parse.unquote_plus(query)
                # logging.debug("Processing token " + str(query))
                for (name, op) in wudb.MyCursor.name_to_operator.items():
                    if op in query:
                        (key, value) = query.split(op, 1)
                        if not name in conditions:
                            conditions[name] = {}
                        # If value is of the form "now(-123)", convert it to a 
                        # time stamp of 123 minutes ago
                        r = re.match(r"now\((-?\d+)\)", value)
                        if r:
                            minutes_ago = int(r.group(1))
                            td = datetime.timedelta(minutes = minutes_ago)
                            value = str(datetime.datetime.now() + td)
                        conditions[name][key] = value
                        break
        logging.debug("send_query(): conditions: %s", conditions)
        if self.db_pool:
            wus = self.db_pool.query(**conditions)
        else:
            wus = wudb.WuAccess(self.dbfilename).query(**conditions)

        body = HtmlGen()

        body.append('<a href="/index.html">Back to index</a>')
        body.append("<p>Query for conditions = " + str(conditions) + "</p>")

        if not wus is None and len(wus) > 0:
            cwd = os.getcwd()
            body.append(str(len(wus)) + " records match.")
            keys = wus[0].keys()
            body.start_table(keys)
            for wu in wus:
                body.wu_row(wu, keys, cwd)
            body.end_table()
        else:
            body.append("No records match.")
        body.finish()
        
        self.send_response(200)
        self.send_header("Content-Type", "text/html")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Content-Length", body.get_len())
        self.end_headers()
        self.send_body(body.__bytes__())


class ServerLauncher(object):
    def __init__(self, address, port, threaded, dbfilename,
                registered_filenames, uploaddir, *, bg = False,
                use_db_pool = True, scriptdir = None, only_registered=False,
                cafile=None, whitelist=None):
        
        self.name = "HTTP server"
        self.logger = logging.getLogger(self.name)
        self.logger.setLevel(logging.NOTSET)
        self.address = address if address else "0.0.0.0"
        self.url_address = address if address else socket.gethostname()
        self.port = port
        self.cafile = cafile
        upload_scriptname = "upload.py"
        # formatter = logging.Formatter(
        #    fmt='%(address_string)s - - [%(asctime)s] %(message)s')
        #self.ch = logging.StreamHandler()
        #self.ch.setFormatter(formatter)
        #self.logger.addHandler(self.ch)
        #self.logger.propagate = False
        
        self.bg = bg
        if use_db_pool:
            self.db_pool = wudb.DbThreadPool(dbfilename, 1)
        else:
            self.db_pool = None
        # Generate a class with parameters stored in class variables
        upload_url_path = "/cgi-bin/upload.py"
        handler_params = {
            "registered_filenames": registered_filenames,
            "logger": self.logger,
            "dbfilename": dbfilename,
            "db_pool": self.db_pool, 
            "uploaddir": uploaddir,
            "cgi_directories" : ['/cgi-bin'],
            "upload_path": upload_url_path,
            "only_registered": only_registered
        }
        MyHandlerWithParams = type("MyHandlerWithParams", (MyHandler, ), handler_params)
        
        # Find the upload.py script
        upload_path = self.findscript(upload_scriptname, scriptdir)
        # Always register the upload script
        if upload_path is None:
            raise IOError("%s script not found" % upload_scriptname)
        self.logger.debug("Found %s at %s" % (upload_scriptname, upload_path))
        registered_filenames[upload_url_path.lstrip('/')] = upload_path

        # Set shell environment variables which the upload.py script needs if
        # spawned as subprocess
        os.environ[upload.DBFILENAMEKEY] = dbfilename
        os.environ[upload.UPLOADDIRKEY] = uploaddir
        if not os.path.isdir(uploaddir):
            os.mkdir(uploaddir)

        # See if we can use HTTPS
        scheme = "http"
        self.cert_sha1 = None
        if self.create_certificate():
            self.cert_sha1 = self.get_certificate_hash()
            if not self.cert_sha1 is None:
                scheme = "https"

        self.url = "%s://%s:%d" % (scheme, self.url_address, self.port)
        
        addr = (self.address, self.port)
        try:
            if threaded and scheme == "http":
                self.logger.info("Using threaded HTTP server")
                self.httpd = ThreadedHTTPServer(addr, MyHandlerWithParams,
                        whitelist=whitelist)
            elif not threaded and scheme == "http":
                self.logger.info("Using non-threaded HTTP server")
                self.httpd = FixedHTTPServer(addr, MyHandlerWithParams,
                        whitelist=whitelist)
            elif threaded and scheme == "https":
                self.logger.info("Using threaded HTTPS server")
                self.httpd = ThreadedHTTPSServer(addr, MyHandlerWithParams,
                        whitelist=whitelist, certfile=self.cafile)
            elif not threaded and scheme == "https":
                self.logger.info("Using non-threaded HTTPS server")
                self.httpd = HTTPSServer(addr, MyHandlerWithParams, 
                        whitelist=whitelist, certfile=self.cafile)
            else:
                assert False
        except socket_error as e:
            if e.errno == errno.EADDRINUSE:
                self.logger.critical("Address %s:%d is already in use (maybe "
                        "another cadofactor running?)", address, port)
                self.logger.critical("You can choose a different port with "
                        "server.port=<integer>.")
                sys.exit(1)
        self.httpd.server_name = self.name

    def get_url(self):
        return self.url

    def get_cert_sha1(self):
        return self.cert_sha1

    def create_certificate(self):
        if self.cafile is None:
            return False
        if os.path.isfile(self.cafile):
            return True
        if not HAVE_SSL:
            self.logger.warn("ssl module not available, won't generate certificate")
            return False
        subj = [
            "C=XY",
            "ST=None",
            "O=None",
            "localityName=None",
            "commonName=%s" % self.url_address,
            "organizationalUnitName=None",
            "emailAddress=None"
        ]
        
        command = ['openssl', 'req', '-new', '-x509', '-batch', '-days', '365',
                   '-nodes', '-subj', '/%s/' % '/'.join(subj), 
                   '-out', self.cafile, '-keyout', self.cafile]
        try:
            output = check_output(command, stderr=STDOUT)
        except (OSError, CalledProcessError) as e:
            self.logger.error("openssl failed: %s", e)
            return False
        self.logger.debug("openssl output: %s", output)
        return True

    def get_certificate_hash(self):
        if self.cafile is None:
            return None
        if not HAVE_SSL:
            self.logger.warn("ssl module not available, won't generate fingerprint")
            return None
        command = ['openssl', 'x509', '-in', self.cafile, '-fingerprint']
        try:
            output = check_output(command)
        except (OSError, CalledProcessError) as e:
            self.logger.error("openssl failed: %s", e)
            return None
        output_text = output.decode("ascii")
        for line in output_text.splitlines():
            if line.startswith("SHA1 Fingerprint="):
                return line.split('=', 1)[1].replace(":", "").lower()
        return None
    
    def serve(self):
        self.logger.info("serving at %s (%s)", self.url, self.httpd.server_address[0])
        certstr = "" if self.cert_sha1 is None else " --certsha1=%s" % self.cert_sha1
        self.logger.info("You can start additional wuclient2.py scripts with "
                         "parameters: --server=%s%s", self.url, certstr)
        self.logger.info("If you want to start additional clients, remember "
                         "to add their hosts to server.whitelist")
        
        if self.bg:
            from threading import Thread
            self.thread = Thread(target=self.httpd.serve_forever,
                                 name=self.name)
            self.thread.daemon = True
            self.thread.start()
        else:
            self.httpd.serve_forever()
    
    def shutdown(self):
        self.logger.info("Shutting down HTTP server")
        self.httpd.shutdown()
        if self.bg:
            self.thread.join()
        if self.db_pool:
            self.db_pool.terminate()
        #self.logger.removeHandler(self.ch)
    
    @staticmethod
    def findscript(scriptname, scriptdir = None):
        # If scriptdir is specified, use that unconditionally
        if not scriptdir is None:
            return scriptdir + os.sep + scriptname
        # Try the CWD
        if os.path.isfile(scriptname):
            return scriptname
        # Try the directory that contains wuserver.py
        if not __file__ is None:
            dirname = os.path.dirname(os.path.realpath(__file__))
            if os.path.isfile(dirname + os.sep + scriptname):
                return dirname + os.sep + scriptname
        # Not found
        return None

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-address", help="Listen address", default="localhost")
    parser.add_argument("-port", help="Listen port", default="8001")
    parser.add_argument("-uploaddir", help="Upload directory", 
                        default="upload/")
    parser.add_argument("-dbfile", help="Database file name", required=True)
    parser.add_argument("-cafile", help="Certificate file name", required=False)
    parser.add_argument("-threaded", help="Use threaded server", 
                        action="store_true", default=False)
    parser.add_argument("-onlyreg", help="Allow access only to registered files", 
                        action="store_true", default=False)
    args = parser.parse_args()

    PORT = int(args.port)
    HTTP = args.address
    dbfilename = args.dbfile
    cafile = args.cafile
    registered_filenames = {}

    logger = logging.getLogger()
    logger.setLevel(logging.NOTSET)

    httpd = ServerLauncher(HTTP, PORT, args.threaded, dbfilename, 
                           registered_filenames, args.uploaddir,
                           only_registered=args.onlyreg, cafile=cafile)
    
    try:
        httpd.serve()
    except KeyboardInterrupt:
        pass
    else:
        raise
    httpd.shutdown()
