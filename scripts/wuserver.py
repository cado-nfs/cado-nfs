#!/usr/bin/env python3

import http.server
import socketserver
import os
import sys
import re
from urllib.parse import unquote_plus
from workunit import Workunit
import datetime
import wudb
import upload

# Get the shell environment variable name in which we should store the path 
# to the upload directory

upload_keywords = ['upload.py']
uploaddir='upload/' # Upload CGI script puts files in this directory
dbfilename='wudb'

class ThreadedHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    """Handle requests in a separate thread."""

debug = 1
def diag(level, text, var = None):
    if debug > level:
        if var is None:
            print (text, file=sys.stderr)
        else:
            print (text + str(var), file=sys.stderr)
        sys.stderr.flush()

class HtmlGen:
    def __init__(self):
        self.body = \
            '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" ' + \
            '"http://www.w3.org/TR/html4/strict.dtd">\n' + \
            '<html>\n' + \
            '<head>\n' + \
            '<title>List of workunits</title>\n' + \
            '</head>\n' + \
            '<body>'

    def __str__(self):
        return self.body + '</body>'

    def append(self, str):
        self.body = self.body + str

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
                    s = s + '<a href="' + path + '">' + f["filename"] + '</a><br>'
                arr.append(s)
            else:
                arr.append(wu[k])
        self.add_table_row(arr)


class MyHandler(http.server.CGIHTTPRequestHandler):
    def send_body(self, body):
        self.wfile.write(bytes(body, "utf-8"))
        self.wfile.flush()

    def do_GET(self):
        """Generates a work unit if request is cgi-bin/getwu, otherwise calls
           parent class' do_GET()"""
        self.cwd = os.getcwd()
        if self.is_cgi():
            if self.is_getwu():
                self.send_WU()
            elif self.is_getstatus():
                self.send_status()
            else:
                self.send_error(404, "GET for CGI scripts allowed only for work unit or status page request")
        else:
            super().do_GET()
        sys.stdout.flush()

    def do_POST(self):
        """Set environment variable telling the upload directory 
           and call CGI handler to run upload CGI script"""
        self.cwd = os.getcwd()
        if self.is_upload():
            os.environ[upload.UPLOADDIRKEY] = uploaddir
            os.environ[upload.DBFILENAMEKEY] = dbfilename
            if False:
                self.send_response(200, "Script output follows")
                upload.do_upload(db, input = self.rfile, output = self.wfile)
            else:
                http.server.CGIHTTPRequestHandler.do_POST(self)
        else:
            self.send_error(404, "POST request allowed only for uploads")
        sys.stdout.flush()

    def is_upload(self):
        """Test whether request is a file upload."""
        splitpath = http.server._url_collapse_path_split(self.path)
        if self.command == 'POST' and self.is_cgi() and splitpath[1] in upload_keywords:
            return True
        return False

    def is_getwu(self):
        """Test whether request is for a new WU."""
        filename=self.cgi_info[1].split("?", 1)[0]
        return self.command == 'GET' and filename in ['getwu']

    def is_getstatus(self):
        """Test whether request is for a a status page."""
        filename=self.cgi_info[1].split("?", 1)[0]
        return self.command == 'GET' and filename in ['status']

    def guess_type(self, path):
        type = super().guess_type(path)
        # Use text/plain for files in upload, unless the type was properly identified
        # FIXME: make path identification more robust
        if type == "application/octet-stream" and path.startswith(self.cwd + '/upload/'):
            return "text/plain"
        return type

    def send_WU(self):
        filename = self.cgi_info[1]
        if not "?" in filename:
            return self.send_error(400, "No query string given")
        (filename, query) = self.cgi_info[1].split("?", 1)
        if query.count("=") != 1 or "?" in query or "&" in query:
            return self.send_error(400, "Bad query string in request")
        (key, clientid) = query.split("=")
        if key != "clientid":
            return self.send_error(400, "No client id specified")
        if not clientid.isalnum():
            return self.send_error(400, "Malformed client id specified")
        
        # wu = wudb.WuAccess(db)
        wu_text = db_pool.assign(clientid)
        if not wu_text:
            return self.send_error(404, "No work available")
        
        self.log_message("Sending work unit " + Workunit(wu_text).get_id() + " to client " + clientid)
        # wu_text = wu.get_wu()
        self.send_response(200)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", len(wu_text))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.send_body(wu_text)

    def send_status(self):
        self.send_query()

    def send_query(self):
        diag(1, "self.cgi_info = ", self.cgi_info)
        filename = self.cgi_info[1]
        if "#" in filename:
            # Get rid of fragment part
            filename = filename.split("#", 1)[0]
        conditions = {}
        if "?" in filename:
            # Parse query part into SELECT conditions
            (filename, query) = filename.split("?", 1)
            print("Query = " + query)
            conditions = {}
            # Now look at individual key=value pairs
            for q in query.split("&"):
                q = unquote_plus(q)
                diag(1, "Processing token ", q)
                for (name, op) in wudb.WuDb.name_to_operator.items():
                    if op in q:
                        (key, value) = q.split(op, 1)
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
        wus = db_pool.query(**conditions)

        body = HtmlGen()

        body.append('<a href="/index.html">Back to index</a>')
        body.append("<p>Query for conditions = " + str(conditions) + "</p>")

        if not wus is None and len(wus) > 0:
            body.append(str(len(wus)) + " records match.")
            keys = wus[0].keys()
            body.start_table(keys)
            for wu in wus:
                body.wu_row(wu, keys, self.cwd)
            body.end_table()
        else:
            body.append("No records match.")
        
        self.send_response(200)
        self.send_header("Content-Type", "text/html")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Content-Length", len(str(body)))
        self.end_headers()
        self.send_body(str(body))

if __name__ == '__main__':
    from sys import argv
    PORT = 8001
    HTTP = "" # address on which to listen
    want_threaded = 1

    if want_threaded:
        ServerClass = ThreadedHTTPServer
    else:
        ServerClass = http.server.HTTPServer

    if argv[1:]:
        PORT = int(argv[1])

    db_pool = wudb.DbThreadPool(dbfilename, 1)

    HandlerClass = MyHandler
    httpd = ServerClass((HTTP, PORT), HandlerClass)
    httpd.server_name = "test"

    print ("serving at " + HTTP + ":" + str(PORT))
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        db_pool.terminate()
    else:
        raise
