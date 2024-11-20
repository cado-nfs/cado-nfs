#!/bin/bash

# This helper script sets up a virtual environment in $PWD/cado-nfs.venv,
# with all the requirements of cado-nfs.

VENV="$PWD/cado-nfs.venv"

if ! [ -d "$VENV" ] ; then
    python3 -m venv $VENV
fi
source $VENV/bin/activate
pip3 install flask requests

echo "To activate this venv, run:"
echo "   source $VENV/bin/activate"
