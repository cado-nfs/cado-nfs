#!/usr/bin/env bash

SCRIPTNAME="$0"

function usage() {
  echo "Usage: $SCRIPTNAME [-s <n>] [-p <path>] [-py <interpreter>] <parameters for wuclient2.py>"
  echo "Mandatory parameters for wuclient2.py: --server=<URL>"
  echo "If the server is configured to use SSL (which is the default), the parameters"
  echo "for wuclient2.py must also include the SSL certificate fingerprint, specified"
  echo "with --certsha1=<CERTSHA1>"
  echo "Optional parameters:"
  echo "-s <n>: Start <n> client scripts on each host. Default: 1"
  echo "-p <path>: Look for wuclient2.py script in directory specified by <path>"
  echo "           on slave nodes. This must be specified if wuclient2.py is not"
  echo "           found in a directory in PATH on the slave nodes."
  echo "-py <interpreter>: Force <interpreter> as the Python interpreter for running"
  echo "                   wuclient2.py on slave nodes. Default: empty."
  echo '-n <n>: Use "numactl --cpunodebind=", cycling through [0, <n>-1] per client'
}


if [ -z "$OAR_NODE_FILE" ]
then
  echo "Error: OAR_NODE_FILE shell environment variable not set." >&2
  echo "$SCRIPTNAME is meant to be run inside an OAR job." >&2
  usage
  exit 1
fi

if [ ! -f "$OAR_NODE_FILE" ]
then
  echo "Error: Variable OAR_NODE_FILE is set, but file $OAR_NODE_FILE does not exist. That's odd." >&2
  usage
  exit 1
fi

# If an -s parameter is given, start that many clients per node

NRCLIENTS=1
SCRIPTPATH=""
unset COMMAND
NUMA_NODE=0
USE_NUMACTL=false
while [ -n "$1" ]
do
  if [ "$1" = "-s" ]
  then
    NRCLIENTS="$2"
    shift 2
  elif [ "$1" = "-p" ]
  then
    SCRIPTPATH="$2"
    shift 2
  elif [ "$1" == "-py" ]
  then
    # Create an array with the Python interpreter, if specified, and later
    # add the path of wuclient2.py. We use an array so that we can expand
    # it to exactly 2 or 1 words, respectively, if a Python interpreter was
    # specified or not, even when any of the paths include whitespace.
    COMMAND[0]="$2"
    shift 2
  elif [ "$1" = "-n" ]
  then
    USE_NUMACTL=true
    NR_CPUS="$2"
    shift 2
  else
    break
  fi
done

# Add trailing slash, if necessary. Note the whitespace before negative offset
# in ${str:offset}
if [[ -n "$SCRIPTPATH" && ! "${SCRIPTPATH: -1}" = "/" ]]
then
  SCRIPTPATH="${SCRIPTPATH}/"
fi
# Add the path of wuclient to COMMAND
COMMAND=("${COMMAND[@]}" "${SCRIPTPATH}wuclient2.py")

if "$USE_NUMACTL"
then
  # Prefix with "numactl" and leave one parameter for --cpunodebind=, to be filled in later
  COMMAND=("numactl" "" "${COMMAND[@]}")
fi

HAVE_SERVER=false
for ARG in "$@"
do
  if [[ "$ARG" =~ "--server" ]]
  then
    HAVE_SERVER=true
  fi
done
if ! "$HAVE_SERVER"
then
  echo "You must specify the URL of the workunit server with --server=<URL>"
  usage
  exit 1
fi


NODES=`cat "$OAR_NODE_FILE" | uniq | tr '\n' ' '`

echo "Starting $NRCLIENTS clients each on $NODES with parameters: $@"

for NODE in $NODES
do 
    for I in `seq 1 "$NRCLIENTS"`
    do
    	if "$USE_NUMACTL"
    	then
    	  # Fill in the --cpunodebind= parameter
    	  COMMAND[1]="--cpunodebind=$NUMA_NODE"
    	  NUMA_NODE="$(( (NUMA_NODE + 1) % NR_CPUS ))"
    	fi
        echo "Running command: oarsh $NODE ${COMMAND[@]} $@"
        oarsh "$NODE" "${COMMAND[@]}" "$@" &
    done
done

# Wait for wuclient jobs to finish.
wait
