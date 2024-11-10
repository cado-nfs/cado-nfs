#!/usr/bin/env python3

# CGI script to handle file uploads

import cgi
import os
import sys
import errno
from tempfile import mkstemp
from shutil import copyfileobj
from cadofactor import wudb

DEBUG = 1

if DEBUG > 1:
    import logging


def diag(level, text, var=None):
    if DEBUG > level:
        if var is None:
            print(text, file=sys.stderr)
        else:
            print("%s%s" % (text, var), file=sys.stderr)
        sys.stderr.flush()


def analyze(level, name, obj):
    """ Dump tons of internal data about an object """
    if DEBUG > level:
        diag(level, "*** Content dump of ", name)
        diag(level, "type(%s): " % name, type(obj))
        diag(level, "dir(%s): " % name, dir(obj))
        diag(level, "%s.__str__() = " % name, obj)
        diag(level, "%s.__repr__() = %r" % (name, obj))
        for name2 in dir(obj):
            diag(level, "%s.%s = " % (name, name2), getattr(obj, name2))


# Global variable in this module so that other Python modules can import
# it and store the path to the upload directory in the shell environment
# variable specified here
UPLOADDIRKEY = "UPLOADDIR"
DBURIKEY = "DBURI"


def do_upload(dburi, uploaddir,
              inputfp=sys.stdin, output=sys.stdout,
              environ=os.environ):
    diag(1, "Command line arguments:", sys.argv)
    diag(2, "Environment:", os.environ)

    try:  # Windows needs stdio set for binary mode.
        # pylint: disable=F0401
        import msvcrt
        # pylint: disable=E1101
        msvcrt.setmode(0, os.O_BINARY)  # stdin  = 0
        msvcrt.setmode(1, os.O_BINARY)  # stdout = 1
    except ImportError:
        pass

    diag(1, "Reading POST data")
    form = cgi.FieldStorage(fp=inputfp, environ=environ)
    diag(1, "Finished reading POST data")
    diag(2, "form = ", form)
    analyze(3, "form", form)

    charset = "utf-8"
    header = "Content-Type: text/plain; charset=" + charset + "\r\n\r\n"

    # A nested FieldStorage instance holds the file
    message = None
    if "WUid" not in form:
        message = 'No "WUid" key found in POST data'
    if "clientid" not in form:
        message = 'No "clientid" key found in POST data'
    elif not os.path.isdir(uploaddir):
        message = 'Script error: %s is not a directory' % uploaddir
    else:
        wuid = form['WUid']
        clientid = form['clientid']
        if 'errorcode' in form:
            errorcode = int(form['errorcode'].value)
            diag(1, "errorcode = ", errorcode)
        else:
            errorcode = None
        if 'failedcommand' in form:
            failedcommand = int(form['failedcommand'].value)
            diag(1, "failedcommand = ", failedcommand)
        else:
            failedcommand = None

        # Test if wuid and clientid was set:
        if not wuid.value:
            message = 'No workunit was specified'
        elif not clientid.value:
            message = 'No client id was specified'

        diag(1, "wuid = ", wuid.value)
        diag(1, "clientid = ", clientid.value)

    if not message:
        filetuples = []
        if 'results' in form:
            fileitems = form['results']
            if isinstance(fileitems, cgi.FieldStorage):
                fileitems = [fileitems]  # Make it iterable
        else:
            fileitems = []
            diag(1, 'No "results" form found')
        analyze(3, "fileitems", fileitems)

        message = ""
        for fileitem in fileitems:
            if not fileitem.file:
                continue
            analyze(3, "f", fileitem)
            diag(1, "Processing file ", fileitem.filename)
            # strip leading path from file name to avoid directory traversal
            # attacks
            basename = os.path.basename(fileitem.filename)
            # Split extension from file name. We need to preserve the
            # file extension so that, e.g., gzipped files can be identified
            # as such
            (basename, suffix) = os.path.splitext(basename)
            # Make a file name which does not exist yet and create the file
            (filedesc, filename) = mkstemp(
                prefix=basename + '.',
                suffix=suffix,
                dir=uploaddir)
            diag(1, "output filename = ", filename)
            filestuple = [fileitem.filename, filename]

            # mkstmp() creates files with mode 0o600 (before umask), and does
            # not allow overriding this with a parameter. We change the mode
            # to 666 & ~umask.

            if os.name != "nt":
                # The os.umask() function gets the old umask *and* sets a new
                # one, so we have to call it twice to avoid changing it :(
                umask = os.umask(0o022)
                os.umask(umask)
                filemode = 0o666 & ~umask
                diag(1, "Setting %s to mode %o" % (filename, filemode))
                os.fchmod(filedesc, filemode)

            filetype = fileitem.headers.get("filetype", None)
            if filetype is not None:
                filestuple.append(filetype)
                diag(1, "filetype = ", filetype)
                command = fileitem.headers.get("command", None)
                if command is not None:
                    filestuple.append(command)
                    diag(1, "command = ", command)
            if False:
                filestuple[1] = os.path.basename(filestuple[1])
            filetuples.append(filestuple)

            # fd is a file descriptor, make a file object from it
            diag(1, "Getting file object for temp file")
            file = os.fdopen(filedesc, "wb")
            diag(1, "Writing data to temp file")
            try:
                copyfileobj(fileitem.file, file)
            except OSError as e:
                if e.errno == errno.ENOSPC:
                    sys.stderr.write("Disk full, deleting %s\n" % filename)
                    file.close()
                    os.unlink(filename)
                raise
            nr_bytes = file.tell()
            diag(1, "Wrote %d bytes" % nr_bytes)
            diag(1, "Closing file")
            file.close()

            # Example output:
            # upload.py: The file "testrun.polyselect.0-5000" for workunit
            # testrun_polyselect_0-5000 was uploaded successfully by client
            # localhost and stored as /localdisk/kruppaal/work/testrun.upload/
            # testrun.polyselect.0-5000.kcudj7, received 84720 bytes.
            message += 'The file "%s" for workunit %s was uploaded ' \
                       'successfully by client %s' \
                       ' and stored as %s, received %d bytes.\n' \
                       % (basename, wuid.value, clientid.value,
                          filename, nr_bytes)
            if errorcode:
                message += 'Error code = %d.\n' % errorcode
        diag(1, "Getting WuAccess object")
        wuar = wudb.WuAccess(dburi)
        diag(1, "Got WuAccess object. Calling .result()")
        try:
            wuar.result(wuid.value, clientid.value, filetuples, errorcode,
                        failedcommand)
        except wudb.StatusUpdateError:
            message = 'Workunit ' + wuid.value + 'was not currently assigned'
        else:
            message = message + 'Workunit ' + wuid.value + ' completed.\n'
        diag(1, "Finished .result()")

    diag(1, sys.argv[0] + ': ', message.rstrip("\n"))

    if output == sys.stdout:
        output.write(header + message)
    else:
        output.write((header + message).encode(charset))


# If this file is run directly by Python, call do_upload()
if __name__ == '__main__':
    if DEBUG > 0:
        import cgitb
        cgitb.enable(display=0, logdir="/tmp/cgitb/")

    for key in (DBURIKEY, UPLOADDIRKEY):
        if key not in os.environ:
            print('Script error: Environment variable %s not set' % key)
            sys.exit(1)

    if DEBUG > 1:
        sys.stderr = open("upload.stderr", "a")
        sys.stderr.write("upload.py: PID = %d\n" % os.getpid())
        logging.basicConfig(level=logging.DEBUG)

    DBURI = os.environ[DBURIKEY]
    UPLOADDIR = os.environ[UPLOADDIRKEY]

    diag(1, "About to call do_upload()")
    do_upload(DBURI, UPLOADDIR)
