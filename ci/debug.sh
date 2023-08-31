#!/usr/bin/env bash


cat <<EOF
# In August 2023, the continuous integration backend of cado-nfs was spun
# off to a separate project (https://gitlab.inria.fr/thome/ci-backend),
# which is now checked out as a submodule of cado-nfs.
#
# In order to use ci/debug.sh, you must now do two things:
#
#  - check out the submodules:
#
#       git submodule update --init --recursive
#
#  - run the ci/ci/debug.sh script, which is the new location of the
#  script you're looking for. (the calling conventions are unchanged)
EOF

exit 1
