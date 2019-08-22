#!/usr/bin/env python3
import os
import sys
import logging

import re
import subprocess
import locale

pathdict=dict()
if not re.search("^/", "@CMAKE_INSTALL_PREFIX@"):
    # We are not in the installed tree, but in the source tree. We need
    # to find the .py libraries.
    pathdict["source"] = os.path.dirname(sys.argv[0])
    pathdict["pylib"] = pathdict["source"] + "/scripts/cadofactor"
    pathdict["data"] = pathdict["source"] + "/parameters"
    # find out where the binaries are. We need to use
    # ./scripts/build_environment.sh for that.
    helper = pathdict["source"] + "/scripts/build_environment.sh"
    pipe = subprocess.Popen([helper, "--show"], stdout=subprocess.PIPE)
    loc = locale.getdefaultlocale()[1]
    if not loc:
        loc="ascii"
    output = pipe.communicate()[0].decode(loc)
    cado_bin_path = [x.split("=",2)[1] for x in output.split("\n") if re.match("^build_tree",x)][0]
    cado_bin_path = re.sub("^\"(.*)\"$", "\\1", cado_bin_path)
    pathdict["lib"] = cado_bin_path
    pathdict["bin"] = cado_bin_path
else:
    pathdict["pylib"] = "@CMAKE_INSTALL_PREFIX@/@LIBSUFFIX@/scripts/cadofactor"
    pathdict["data"] ="@CMAKE_INSTALL_PREFIX@/@DATASUFFIX@"
    # binaries are installed in subdirectories of $LIBDIR.
    pathdict["lib"] ="@CMAKE_INSTALL_PREFIX@/@LIBSUFFIX@"
    pathdict["bin"] ="@CMAKE_INSTALL_PREFIX@/@BINSUFFIX@"

# note that even though we do have cado-nfs.py and cado-nfs-client.py in
# the build tree, we make *NO PROMISE* as to whether calling these
# scripts works. Only calling either within the source tree or within the
# installed tree is expected to work.
sys.path.append(pathdict["pylib"])


import cadotask
import cadologger
import toplevel
import itertools
import wudb
from cadocommand import shellquote

if __name__ == '__main__':
    # Parse command line arguments
    
    # Some command-line arguments are really parsed only here, while some
    # others are relevant to the whole hierarchy of cado-nfs programs.
    # The (hairy) logic which is used to form the definitive list of
    # parameters (in the cadoparams sense) from what we got here on the
    # command line is grouped in the Cado_NFS_toplevel class, down in
    # scripts/cadofactor/toplevel.py

    toplevel_params = toplevel.Cado_NFS_toplevel()
    for key, value in pathdict.items():
        toplevel_params.setpath(key, value)
    logger = toplevel_params.logger
    parameters, db = toplevel_params.get_cooked_parameters()

    # We have to remove the credential info from the database which gets
    # stored in the parameter snapshot. We must also make sure that we
    # store this info always at the root of the param tree, or we may end
    # up do bizarre things with duplicated keys if we resume from a
    # parameter snapshot file *and* we have something on the command
    # line.
    if parameters.get_or_set_default("database", None):
        parameters.replace("database", db.uri_without_credentials)
        # do this so that the parameter does not appear unused.
        parameters.get_or_set_default("database")


    # well, this *must* exist, right ?
    name = parameters.get_or_set_default("tasks.name")
    wdir = parameters.get_or_set_default("tasks.workdir")
    
    # Add a logger to capture the command lines of programs we run
    cmdfilename = os.path.join(wdir, name + ".cmd")
    logger.addHandler(cadologger.CmdFileHandler(cmdfilename))
    
    # Add a logger to write debugging information to a log file
    filelvl = getattr(cadologger, toplevel_params.args.filelog.upper())
    logfilename = os.path.join(wdir, name + ".log")
    filehandler = cadologger.FileHandler(filename = logfilename, lvl = filelvl)
    logger.addHandler(filehandler)

    
    cmdline=" ".join([shellquote(arg, idx == 0) for idx, arg in enumerate(sys.argv)])
    cmdline=cmdline.replace(db.uri, db.uri_without_credentials)
    logger.info("Command line parameters: %s", cmdline)


    logger.debug("Root parameter dictionary:\n%s", parameters)


    # Write a snapshot of the parameters to a file
    for counter in itertools.count():
        snapshot_basename = name + ".parameters_snapshot.%d" % counter
        snapshot_filename = os.path.join(wdir, snapshot_basename)
        if not os.path.isfile(snapshot_filename):
            break
    with open(snapshot_filename, "w") as snapshot_file:
        logger.debug("Writing parameter snapshot to %s", snapshot_filename)
        snapshot_file.write(str(parameters))
        snapshot_file.write("\n")

    logger.info("If this computation gets interrupted, it can be resumed with %s %s", sys.argv[0], snapshot_filename)

    factorjob = cadotask.CompleteFactorization(db=db,
                                               parameters = parameters,
                                               path_prefix = [])

    if toplevel_params.args.verboseparam:
        logger.info("Summary of all recognized parameters\n" +
                factorjob.parameter_help)

    factors = factorjob.run()
    
    dlp_param = parameters.myparams({"dlp": False,}, "")
    dlp = dlp_param["dlp"]
    checkdlp_param = parameters.myparams({"checkdlp": True ,}, "")
    checkdlp = checkdlp_param["checkdlp"]
    target_param = parameters.myparams({"target": "",}, "")
    target = target_param["target"]
    if factors is None:
        toplevel_params.purge_temp_files(nopurge=True)
        sys.exit("Error occurred, terminating")
    else:
        toplevel_params.purge_temp_files()

    if not dlp:
        print(" ".join(factors))
    else:
        if checkdlp:
            p = int(factors[0])
            ell = int(factors[1])
            log2 = int(factors[2])
            log3 = int(factors[3])
            logger.info("Checking that log(2) and log(3) are consistent...")
            logger.info("  p = " + str(p))
            logger.info("  ell = " + str(ell))
            logger.info("  log2 = " + str(log2))
            logger.info("  log3 = " + str(log3))
            assert (p-1) % ell == 0
            assert pow(3, log2*((p-1) // ell), p) == pow(2, log3*((p-1) // ell), p)
            if target != 0:
                logtarget = int(factors[4])
                logger.info("Also check log(target) vs log(2) ...")
                assert pow(int(target), log2*((p-1) // ell), p) == pow(2, logtarget*((p-1) // ell), p)
        else:
            logger.info("No check was performed. Logarithms of the factor base elements are in %s" % factorjob.request_map[cadotask.Request.GET_DLOG_FILENAME]())
        if target != 0:
            logtarget = int(factors[4])
            logger.info("target = " + str(target))
            logger.info("log(target) = " + str(logtarget))
            print(str(logtarget))
        logger.info("If you want to compute a new target, run %s %s target=<target>", sys.argv[0], snapshot_filename)

