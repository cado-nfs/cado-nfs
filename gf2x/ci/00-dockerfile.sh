#!/bin/sh

# We're /bin/sh, not bash.
#
# Our output must be a Dockerfile

HUSH_STDOUT=1

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

FROM=alpine:latest

if [ "$clang" ] ; then
    # difficult. fedora:rawhide, debian:testing are ok and seem to bring
    # us the latest clang. But they're both quite heavy!
    # FROM=debian:testing
    # This is just as heavy (if not worse), but does bring the benefit
    # that we have the latest clang, as is the case with gcc.
    FROM=silkeh/clang:latest
fi
if [ "$gcc" ] && ! [ "$coverage" ] ; then
    # coverage tests need gcc, but also need recent gcov. We can get by
    # with the gcc on alpine linux.
    FROM=gcc:latest
    # well, thinking of it... alpine linux, as of today, is on gcc
    # preversions...
fi

if [ "$icc" ] ; then
    FROM=intel/oneapi-hpckit:latest
fi

# it's important that these come after the per-compiler selection of
# "ideal" images, of course.
#
# note that we're not matching the "on " in "on debian9 system", so that
# an argument "debian9 system with clang" should just work, which sounds
# a bit more natural than keeping "checks on " in the argument.
case "$CI_JOB_NAME" in
    *"alpine system"*) FROM=alpine:latest;;
    *"alpine-edge system"*) FROM=alpine:edge;;
    *"debian system"*) FROM=debian;;
    *"opensuse system"*) FROM=opensuse/leap;;
    *"debian8 system"*) FROM=debian:8;;
    *"debian9 system"*) FROM=debian:9;;
    *"debian10 system"*) FROM=debian:10;;
    *"debian11 system"*) FROM=debian:11;;
    *"debian-testing system"*) FROM=debian:testing;;
    *"debian-unstable system"*) FROM=debian:unstable;;
    *"fedora25 system"*) FROM=fedora:25;;
    *"fedora26 system"*) FROM=fedora:26;;
    *"fedora27 system"*) FROM=fedora:27;;
    *"fedora28 system"*) FROM=fedora:28;;
    *"fedora29 system"*) FROM=fedora:29;;
    *"fedora30 system"*) FROM=fedora:30;;
    *"fedora31 system"*) FROM=fedora:31;;
    *"fedora32 system"*) FROM=fedora:32;;
    *"fedora33 system"*) FROM=fedora:33;;
    *"fedora34 system"*) FROM=fedora:34;;
    *"fedora-rawhide system"*) FROM=fedora:rawhide;;
esac

cat <<EOF
FROM $FROM
EOF

if [ "$DOCKER_SCRIPT" ] ; then
    echo "ENV DOCKER_SCRIPT=1"
fi

if [ "$CI_JOB_NAME" ] ; then
    echo "ENV CI_JOB_NAME=\"$CI_JOB_NAME\""
fi

if [ "$1" = debug ] ; then
cat <<EOF
COPY ./ /tmp/ci/
RUN bash -x /tmp/ci/00-prepare-docker.sh
EOF
else
cat <<EOF
COPY ./ /tmp/ci/
RUN /tmp/ci/00-prepare-docker.sh
EOF
fi

if [ "$icc" ] ; then
    cat <<EOF
RUN ln -s $ONEAPI_ROOT/setvars.sh /etc/profile.d/90-intel-setvars.sh
EOF
fi
