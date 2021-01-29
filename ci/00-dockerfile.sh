#!/bin/sh

# We're /bin/sh, not bash.
#
# Our output must be a Dockerfile

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

FROM=alpine

if [ "$clang" ] ; then
    # difficult. fedora:rawhide, debian:testing are ok and seem to bring
    # us the latest clang. But they're both quite heavy!
    FROM=debian:testing
fi
if [ "$gcc" ] && ! [ "$coverage" ] ; then
    # coverage tests need gcc, but also need recent gcov. We can get by
    # with the gcc on alpine linux.
    FROM=gcc
    # well, thinking of it... alpine linux, as of today, is on gcc
    # preversions...
fi

if [ "$icc" ] ; then
    FROM=intel/oneapi-hpckit
fi

# it's important that these come after the per-compiler selection of
# "ideal" images, of course.
case "$CI_BUILD_NAME" in
    *"on alpine system"*) FROM=alpine;;
    *"on alpine-edge system"*) FROM=alpine:edge;;
    *"on debian8 system"*) FROM=debian:8;;
    *"on debian9 system"*) FROM=debian:9;;
    *"on debian10 system"*) FROM=debian:10;;
    *"on debian11 system"*) FROM=debian:11;;
    *"on debian-testing system"*) FROM=debian:testing;;
    *"on debian-unstable system"*) FROM=debian:unstable;;
    *"on fedora25 system"*) FROM=fedora:25;;
    *"on fedora26 system"*) FROM=fedora:26;;
    *"on fedora27 system"*) FROM=fedora:27;;
    *"on fedora28 system"*) FROM=fedora:28;;
    *"on fedora29 system"*) FROM=fedora:29;;
    *"on fedora30 system"*) FROM=fedora:30;;
    *"on fedora31 system"*) FROM=fedora:31;;
    *"on fedora32 system"*) FROM=fedora:32;;
    *"on fedora33 system"*) FROM=fedora:33;;
    *"on fedora34 system"*) FROM=fedora:34;;
    *"on fedora-rawhide system"*) FROM=fedora:rawhide;;
esac

cat <<EOF
FROM $FROM
EOF

if [ "$DOCKER_SCRIPT" ] ; then
    echo "ENV DOCKER_SCRIPT=1"
fi

if [ "$CI_BUILD_NAME" ] ; then
    echo "ENV CI_BUILD_NAME=\"$CI_BUILD_NAME\""
fi

cat <<EOF
COPY 0??-*.sh ??-*.sh /tmp/
RUN /tmp/00-prepare-docker.sh
EOF

if [ "$icc" ] ; then
    cat <<EOF
RUN ln -s $ONEAPI_ROOT/setvars.sh /etc/profile.d/90-intel-setvars.sh
EOF
fi
