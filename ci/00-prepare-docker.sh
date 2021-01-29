#!/bin/sh

# We're /bin/sh, not bash.
#
# Our output must be a Dockerfile

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

if ! [ "$DOCKER_SCRIPT" ] ; then
    echo "Enter CI script for $CI_PROJECT_NAMESPACE/$CI_PROJECT_NAME, stage $CI_JOB_STAGE ; $CI_BUILD_NAME"
fi

enter_section preparation "System preparation (docker)"
# id -a
# use this to dump environment variables
# export
echo "en_US.UTF-8 UTF-8" > /etc/locale.gen
leave_section

enter_section install_packages "Installing required packages"

debian_packages="$debian_packages     bc"
debian_packages="$debian_packages     locales"
debian_packages="$debian_packages     cmake"
debian_packages="$debian_packages     libhwloc-dev"
debian_packages="$debian_packages     libgmp-dev"
debian_packages="$debian_packages     python3"

fedora_packages="$fedora_packages     bc"
fedora_packages="$fedora_packages     cmake"
fedora_packages="$fedora_packages     hwloc-devel"
fedora_packages="$fedora_packages     gmp-devel"
fedora_packages="$fedora_packages     hostname"
fedora_packages="$fedora_packages     python"

alpine_packages="$alpine_packages     bc"
alpine_packages="$alpine_packages     cmake"
alpine_packages="$alpine_packages     hwloc-dev"
alpine_packages="$alpine_packages     gmp-dev"
alpine_packages="$alpine_packages     make"
alpine_packages="$alpine_packages     bash"
alpine_packages="$alpine_packages     perl"
alpine_packages="$alpine_packages     python3"

while [ $# -gt 0 ] ; do
    case "$1" in
        coverage|clang|gcc|debug|icc) eval "$1=1";;
        *) echo "$1 -> ???" >&2 ; exit 1;;
    esac
    shift
done

# These variables are set in ci/001-environment.sh
if [ "$coverage" ] ; then
    debian_packages="$debian_packages     gcovr"
    fedora_packages="$fedora_packages     gcovr"
    alpine_packages="$alpine_packages     gcovr"
fi

if [ "$gcc" ] ; then
    debian_packages="$debian_packages     g++"
    fedora_packages="$fedora_packages     g++"
    alpine_packages="$alpine_packages     g++"
fi

if [ "$clang" ] ; then
    debian_packages="$debian_packages     clang"
    fedora_packages="$fedora_packages     clang"
    alpine_packages="$alpine_packages     clang"
fi

if [ "$DOCKER_SCRIPT" ] ; then
    debian_packages="$debian_packages sudo git vim gdb"
    fedora_packages="$fedora_packages sudo git vim gdb"
    alpine_packages="$alpine_packages sudo git vim gdb"
fi

if [ -f /etc/debian_version ] ; then
    DEBIAN_FRONTEND=noninteractive apt-get -y update
    DEBIAN_FRONTEND=noninteractive apt-get -y install $debian_packages
elif [ -f /etc/fedora-release ] ; then
    dnf -y install $fedora_packages
elif [ -f /etc/alpine-release ] ; then
    # hwloc-dev still in alpine testing.
    cat >> /etc/apk/repositories <<EOF
http://dl-cdn.alpinelinux.org/alpine/edge/testing
EOF
    apk update
    apk add $alpine_packages
fi

leave_section

