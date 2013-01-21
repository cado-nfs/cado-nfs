#!/usr/bin/env python3

import sys
import os
import stat
import argparse
import shutil
import time
import urllib.request
import subprocess
import hashlib
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import email.encoders
import email.generator
from string import Template
from io import StringIO,BytesIO
from workunit import Workunit

class FixedBytesGenerator(email.generator.BytesGenerator):
    def _handle_bytes(self, msg):
        payload = msg.get_payload()
        if payload is None:
            return
        if isinstance(payload, bytes):
            # Payload is bytes, output is bytes - just write them
            self._fp.write(payload)
        else:
            # Payload is neither bytes not string - this can't be right
            raise TypeError('bytes payload expected: %s' % type(payload))
    _writeBody = _handle_bytes


# Settings which we require on the command line (no defaults)
REQUIRED_SETTINGS = {"CLIENTID" : ("", "Unique ID for this client"), 
                     "DLDIR" : ("", "Directory for downloading files"), 
                     "SERVER" : ("", "Base URL for WU server"), 
                     "WORKDIR" : ("", "Directory for result files")}
# Optional settings with defaults, overrideable on command line, and a help text
OPTIONAL_SETTINGS = {"WU_FILENAME" : ("WU", "Filename under which to store WU files"), 
                     "GETWUPATH" : ("/cgi-bin/getwu", "Path segment of URL for requesting WUs from server"), 
                     "POSTRESULTPATH" : ("/cgi-bin/upload.py", "Path segment of URL for reporting results to server"), 
                     "DEBUG" : ("0", "Debugging verbosity"),
                     "ARCH" : ("", "Architecture string for this client"),
                     "DOWNLOADRETRY" : ("300", "Time to wait before download retries"),
                     "NICENESS" : ("0", "Run subprocesses under this niceness")}
# Merge the two, removing help string
SETTINGS = dict([(a,b) for (a,(b,c)) in list(REQUIRED_SETTINGS.items()) + \
                                        list(OPTIONAL_SETTINGS.items())])

def log(level, text, var = None):
    if int(SETTINGS["DEBUG"]) >= level:
        if var == None:
            print (text)
        else:
            print (text + str(var))

def get_file(urlpath, dlpath = None, options = None):
    # print('get_file("' + urlpath + '", "' + dlpath + '")')
    if dlpath == None:
        dlpath = SETTINGS["DLDIR"] + "/" + os.basename(urlpath) # FIXME: should be url base name
    url = SETTINGS["SERVER"] + "/" + urlpath
    if options:
        url = url + "?" + options
    log (0, "Downloading " + url + " to " + dlpath);
    request = None
    while request == None:
        try:
            request = urllib.request.urlopen(url)
        except urllib.error.HTTPError as e:
            log (0, str(e))
            return False
        except (NameError, urllib.error.URLError) as e:
            request = None
            wait = float(SETTINGS["DOWNLOADRETRY"])
            log(0, "Download of " + urlpath + " failed, " + str(e))
            log(0, "Waiting " + str(wait) + " seconds before retrying")
            time.sleep(wait)
    file = open(dlpath, "wb")
    shutil.copyfileobj (request, file)
    file.close()
    request.close()
    return True

def do_checksum(filename, checksum = None):
    """ Computes the SHA1 checksum for a file. If checksum is None, returns 
        the computed checksum. If checksum is not None, return whether the
        computed SHA1 sum and checksum agree """
    blocksize = 65536
    m = hashlib.sha1()
    file = open(filename, "rb")
    data = file.read(blocksize)
    while data:
        m.update(data)
        data = file.read(blocksize)
    file.close()
    filesum = m.hexdigest()
    if checksum is None:
        return filesum
    else:
        return filesum.lower() == checksum.lower()

def get_missing_file(urlpath, filename, checksum = None, options = None):
    # print('get_missing_file("' + urlpath + '", "' + filename + '", ' + str(checksum) + ')')
    if os.path.isfile(filename):
        log (0, filename + " already exists, not downloading")
        if checksum is None:
            return True
        filesum = do_checksum(filename)
        if filesum.lower() == checksum.lower():
            return True
        log (0, "Existing file " + filename + " has wrong checksum " + filesum + 
             ", workunit specified " + checksum +". Deleting file.")
        os.remove(filename)
    
    # If checksum is wrong and does not change during two downloads, exit with 
    # failue, as apparently the file on the server and checksum in workunit do 
    # not agree
    last_filesum = None
    while True:
        if not get_file(urlpath, filename, options):
            return False
        if checksum is None:
            return True
        filesum = do_checksum(filename)
        if filesum.lower() == checksum.lower():
            return True
        if not last_filesum is None and filesum == last_filesum:
            log (0, "Downloaded file " + filename + " has same wrong checksum " 
                 + filesum + " again. Exiting.")
            return False
        log (0, "Downloaded file " + filename + " has wrong checksum " + 
             filesum + ", workunit specified " + checksum + ". Deleting file.")
        os.remove(filename)
        last_filesum = filesum

class Workunit_Processor(Workunit):
    
    def __init__(self, filepath, debug = 0):
        self.exitcode = 0 # will get set if any command exits with code != 0
        self.failedcommand = None
        self.clientid = SETTINGS["CLIENTID"]
        self.debug = debug # Controls debugging output

        log (1, "Parsing workunit from file " + filepath)
        wu_file = open(filepath)
        wu_text = wu_file.read()
        wu_file.close()
        self.wu = Workunit(wu_text)
        self.WUid = self.wu.get_id()
        log (1, " done, workunit ID is " + self.WUid)
        self.stdout = []
        self.stderr = []

    def  __str__(self):
        return "Processor for Workunit:\n" + self.wu.__str__()

    def get_files(self):
        for (filename, checksum) in self.wu.data["FILE"] + self.wu.data["EXECFILE"]:
            archname = Template(filename).safe_substitute({"ARCH": SETTINGS["ARCH"]})
            dlname = Template(filename).safe_substitute({"ARCH": ""})
            if not get_missing_file (archname, SETTINGS["DLDIR"] + '/' + dlname, checksum):
                return False
        for (filename, checksum) in self.wu.data["EXECFILE"]:
            dlname = Template(filename).safe_substitute({"ARCH": ""})
            path = SETTINGS["DLDIR"] + '/' + dlname
            mode = os.stat(path).st_mode
            if mode & stat.S_IXUSR == 0:
                log (0, "Setting executable flag for " + path)
                os.chmod(path, mode | stat.S_IXUSR)
        return True

    @staticmethod
    def renice():
        os.nice(int(SETTINGS["NICENESS"]))
    
    def run_commands(self):
        for (counter, command) in enumerate(self.wu.data["COMMAND"]):
            command = Template(command).safe_substitute(SETTINGS)
            log (0, "Running command for workunit " + self.wu.get_id() + ": " + command)

            # If niceness command line parameter was set, call self.renice() in
            # child process, before executing command
            if int(SETTINGS["NICENESS"]) > 0:
                renice_func = self.renice
            else:
                renice_func = None

            # Run the command
            child = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, 
                stderr = subprocess.PIPE, preexec_fn = renice_func)
            # Wait for command to finish executing, capturing stdout and stderr 
            # in output tuple
            (child_stdout, child_stderr) = child.communicate()

            if len(child_stdout) > 0:
                self.stdout.append(child_stdout)
            else:
                self.stdout.append(None)
            if len(child_stderr) > 0:
                self.stderr.append(child_stderr)
            else:
                self.stderr.append(None)

            if child.returncode != 0:
                log (0, "Command exited with exit code " + str(child.returncode)) 
                self.failedcommand = counter
                self.exitcode = child.returncode
                return False
            else:
                log (1, "Command exited successfully")
        return True

    def upload_result(self):
        # Build a multi-part MIME document containing the WU id and result file
        postdata = MIMEMultipart()
        for key in ("WUid", "clientid", "exitcode", "failedcommand"):
            if not getattr(self, key, None) is None:
                attachment = MIMEText(str(getattr(self, key)))
                attachment.add_header('Content-Disposition', 'form-data', name=key)
                postdata.attach(attachment)
        if "RESULT" in self.wu.data:
            for f in self.wu.data["RESULT"]:
                filepath = SETTINGS["WORKDIR"] + "/" + f
                log (1, "Adding result file " + filepath + " to upload")
                file = open(filepath, "rb")
                filedata = file.read()
                file.close()
                # HTTP does not use a Content-Transfer-Encoding, so use noop encoder
                # This fails, probably related to http://bugs.python.org/issue4768
                result = MIMEApplication(filedata, _encoder=email.encoders.encode_noop)
                result.add_header('Content-Disposition', 'form-data', name="results", 
                                  filename=f)
                postdata.attach(result)
        for (name, arr) in (("stdout", self.stdout), ("stderr", self.stderr)):
            for (counter, f) in enumerate(arr):
                if not f is None:
                    log (1, "Adding " + name + " for command " + str(counter) + " to upload")
                    result = MIMEApplication(f, _encoder=email.encoders.encode_noop)
                    result.add_header('Content-Disposition', 'form-data', name="results", 
                                      filename=name + str(counter))
                    postdata.attach(result)
        
        if self.debug >= 2:
            print("Headers of postdata as a dictionary:")
            print(dict(postdata.items()))
        # Ugly hack: overwrite method for writing headers to suppress them
        # We pass the MIME headers to the request below via the headers= argument, 
        # and don't want them to occur again as part of the POST data
        if False:
            postdata2 = postdata.as_string(unixfrom=False) + "\n"
            if self.debug >= 2:
                print("Postdata as a string:")
                print(postdata2)
            postdata3 = bytes(postdata2, encoding="utf-8")
        else:
            fp = BytesIO()
            g = FixedBytesGenerator(fp)
            g.flatten(postdata, unixfrom=False)
            postdata3 = fp.getvalue() + b"\n"
        if self.debug >= 2:
            print("Postdata as a bytes array:")
            print(postdata3)
        url = SETTINGS["SERVER"] + "/" + SETTINGS["POSTRESULTPATH"]
        request = urllib.request.Request(url, data=postdata3, headers=dict(postdata.items()))
        conn = None;
        while conn is None:
            try:
                conn = urllib.request.urlopen(request)
            except (urllib.error.URLError) as e:
                conn = None
                wait = float(SETTINGS["DOWNLOADRETRY"])
                log(0, "Upload of result failed, " + str(e))
                log(0, "Waiting " + str(wait) + " seconds before retrying")
                time.sleep(wait)
        log (1, "Server response:")
        for line in conn:
            log(1, line)
        conn.close()
        return True

    def result_exists(self):
        if "RESULT" in self.wu.data:
            for f in self.wu.data["RESULT"]:
                filepath = SETTINGS["WORKDIR"] + "/" + f
                if not os.path.isfile(filepath):
                    log (0, "Result file " + filepath + " does not exist")
                    return False
                log (0, "Result file " + filepath + " already exists")
        log (0, "All result files already exist")
        return True

    def cleanup(self):
        log (0, "Cleaning up for workunit " + self.wu.get_id())
        if "RESULT" in self.wu.data:
            for f in self.wu.data["RESULT"]:
                filepath = SETTINGS["WORKDIR"] + "/" + f
                log (0, "Removing result file " + filepath)
                os.remove(filepath)

    def process(self):
        # If all output files exist, send them, return WU as done
        # Otherwise, run commands in WU. If no error and all output 
        #   files exist, send them, return WU as done
        # print(str(wu))
        if not self.get_files():
            return False
        if not self.result_exists():
            if not self.run_commands():
                return False
        if not self.upload_result():
            return False
        self.cleanup()
        return True

def do_work():
    wu_filename = SETTINGS["DLDIR"] + "/" + SETTINGS["WU_FILENAME"]
    if not get_missing_file(SETTINGS["GETWUPATH"], wu_filename, 
                            options="clientid=" + SETTINGS["CLIENTID"]):
        return False
    wu = Workunit_Processor(wu_filename, int(SETTINGS["DEBUG"]))
    if not wu.process():
        return False
    log (0, "Removing workunit file " + wu_filename)
    os.remove(wu_filename)
    return True

if __name__ == '__main__':
    # Create command line parser from the keys in SETTINGS
    parser = argparse.ArgumentParser()
    for arg in REQUIRED_SETTINGS.keys():
        parser.add_argument('--' + arg.lower(), required = True,
        help=REQUIRED_SETTINGS[arg][1])
    for arg in OPTIONAL_SETTINGS.keys():
        parser.add_argument('--' + arg.lower(), required = False, 
            default=OPTIONAL_SETTINGS[arg][0], 
            help=OPTIONAL_SETTINGS[arg][1] + " (default: " + OPTIONAL_SETTINGS[arg][0] + ")")
    # Parse command line, store as dictionary
    args = vars(parser.parse_args())
    # Copy values to SETTINGS
    for arg in SETTINGS.keys():
        if arg.lower() in args:
            SETTINGS[arg] = args[arg.lower()]

    # print (str(SETTINGS))

    while do_work():
        pass

