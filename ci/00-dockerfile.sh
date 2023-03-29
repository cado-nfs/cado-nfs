#!/bin/sh

# We're /bin/sh, not bash.
#
# Our output must be a Dockerfile

HUSH_STDOUT=1

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/005-build-environment.sh"

FROM=alpine:latest

if [ "$clang" ] ; then
    # difficult. fedora:rawhide, debian:testing are ok and seem to bring
    # us the latest clang. But they're both quite heavy!
    # FROM=debian:testing
    # This is just as heavy (if not worse), but does bring the benefit
    # that we have the latest clang, as is the case with gcc.
    FROM=silkeh/clang:latest
    case "$clang" in
        dev) FROM=silkeh/clang:dev;;
        [0-9][0-9]*) FROM=silkeh/clang:$clang;;
    esac
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
case "$BUILD_NAME" in
    *"alpine system"*) FROM=alpine:latest;;
    *"alpine-edge system"*) FROM=alpine:edge;;
    *"debian system"*) FROM=debian:latest;;
    *"opensuse system"*) FROM=opensuse/leap;;
    *"fedora system"*) FROM=fedora:latest;;
    *"centos system"*) FROM=quay.io/centos/centos:stream;;
    *"centos9 system"*) FROM=quay.io/centos/centos:stream9;;
    *"debian8 system"*) FROM=debian:8;;
    *"debian9 system"*) FROM=debian:9;;
    *"debian10 system"*) FROM=debian:10;;
    *"debian11 system"*) FROM=debian:11;;
    *"debian-testing system"*) FROM=debian:testing;;
    *"debian-unstable system"*) FROM=debian:unstable;;
    *"ubuntu system"*) FROM=ubuntu:latest;;
    *"ubuntu rolling system"*) FROM=ubuntu:rolling;;
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
    *"containers"*) FROM=docker;;
esac

if [ "$1" = "--get-container" ] ; then
    echo $FROM
    exit 0
fi

cat <<EOF
FROM $FROM
EOF

if [ "$DOCKER_SCRIPT" ] ; then
    echo "ENV DOCKER_SCRIPT=1"
fi

if [ "$BUILD_NAME" ] ; then
    echo "ENV BUILD_NAME=\"$BUILD_NAME\""
fi

cat <<EOF
COPY ./ /tmp/ci/
RUN /tmp/ci/00-prepare-docker.sh
EOF

if [ "$icc" ] ; then
    cat <<EOF
RUN ln -s $ONEAPI_ROOT/setvars.sh /etc/profile.d/90-intel-setvars.sh
EOF
fi
