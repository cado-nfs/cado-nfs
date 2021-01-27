#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

if ! [ "$DOCKER_SCRIPT" ] ; then
    echo "Enter CI script for $CI_PROJECT_NAMESPACE/$CI_PROJECT_NAME, stage $CI_JOB_STAGE ; $CI_BUILD_NAME"
fi

enter_section preparation "System preparation (docker)"
id -a
# use this to dump environment variables
# export
echo "en_US.UTF-8 UTF-8" > /etc/locale.gen
leave_section

enter_section install_packages "Installing required packages"

debian_packages=(
    bc
    locales
    cmake
    libhwloc-dev
    libgmp-dev
)

fedora_packages=(
    bc
    cmake
    hwloc-devel
    gmp-devel
    hostname
)

while [ $# -gt 0 ] ; do
    case "$1" in
        coverage|clang|gcc|debug|icc) eval "$1=1";;
        *) echo "$1 -> ???" >&2 ; exit 1;;
    esac
    shift
done

if [[ $CI_BUILD_NAME =~ "coverage tests" ]] || [ "$coverage" ] ; then
    debian_packages+=(gcovr)
    fedora_packages+=(gcovr)
fi

if [[ $CI_BUILD_NAME =~ "with gcc" ]] || [ "$gcc" ] ; then
    debian_packages+=(g++)
    fedora_packages+=(g++)
fi

if [[ $CI_BUILD_NAME =~ "with clang" ]] || [ "$clang" ] ; then
    debian_packages+=(clang)
    fedora_packages+=(clang)
fi

if [ "$DOCKER_SCRIPT" ] ; then
    debian_packages+=(sudo python3 git vim gdb g++)
    fedora_packages+=(sudo python git vim gdb gcc g++ make)
fi

if [ -f /etc/debian_version ] ; then
    DEBIAN_FRONTEND=noninteractive apt-get -y update
    DEBIAN_FRONTEND=noninteractive apt-get -y install "${debian_packages[@]}"
elif [ -f /etc/fedora-release ] ; then
    dnf -y install "${fedora_packages[@]}"
fi

leave_section

