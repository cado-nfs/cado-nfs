#!/bin/sh

# We're /bin/sh, not bash.
#
# Our output must be a Dockerfile

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

if ! [ "$DOCKER_SCRIPT" ] ; then
    echo "Enter CI script for $CI_PROJECT_NAMESPACE/$CI_PROJECT_NAME, stage $CI_JOB_STAGE ; $CI_JOB_NAME"
fi

enter_section preparation "System preparation (${RUNTIME_TYPE:-docker})"
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
# debian has gmp-ecm. Let's use it in our debian tests, that will give
# some extra coverage.
debian_packages="$debian_packages     libecm-dev"
# is full perl really needed ?
# debian_packages="$debian_packages     perl"
debian_packages="$debian_packages     python3"
# we may consider using a system libfmt at some point.
# debian_packages="$debian_packages     libfmt-dev"

opensuse_packages="$opensuse_packages     bc"
opensuse_packages="$opensuse_packages     which"  # is type -p more portable?
opensuse_packages="$opensuse_packages     hostname"
opensuse_packages="$opensuse_packages     cmake"
opensuse_packages="$opensuse_packages     hwloc-devel"
opensuse_packages="$opensuse_packages     gmp-devel"
opensuse_packages="$opensuse_packages     python3"
opensuse_packages="$opensuse_packages     gzip"

fedora_packages="$fedora_packages     util-linux-core"
fedora_packages="$fedora_packages     bc"
fedora_packages="$fedora_packages     cmake"
fedora_packages="$fedora_packages     hwloc-devel"
fedora_packages="$fedora_packages     gmp-devel"
fedora_packages="$fedora_packages     hostname"
# is full perl really needed ? seems that perl-interpreter and the auto
# dependencies that we pull already pull what we need.
fedora_packages="$fedora_packages     perl-interpreter"
fedora_packages="$fedora_packages     python"
fedora_packages="$fedora_packages     findutils diffutils"

centos_packages="$centos_packages     bc"
centos_packages="$centos_packages     cmake"
centos_packages="$centos_packages     hwloc"
if is_centos_above_8 ; then
    # alright -- centos stream 8 does not have hwloc-devel. It's just
    # weird...
    centos_packages="$centos_packages     hwloc-devel"
fi
centos_packages="$centos_packages     gmp-devel"
centos_packages="$centos_packages     hostname"
centos_packages="$centos_packages     perl-interpreter"
centos_packages="$centos_packages     python3"
centos_packages="$centos_packages     findutils diffutils"

alpine_packages="$alpine_packages     bc"
alpine_packages="$alpine_packages     cmake"
# alpine_packages="$alpine_packages     hwloc-dev"
alpine_packages="$alpine_packages     gmp-dev"
alpine_packages="$alpine_packages     make"
alpine_packages="$alpine_packages     bash"
alpine_packages="$alpine_packages     perl"
alpine_packages="$alpine_packages     python3"
alpine_packages="$alpine_packages     gzip"

freebsd_packages="$freebsd_packages     cmake"
# See #30036. We NEVER want to include hwloc under freebsd.
# 	hwloc: 1.11.13_1 seems to be the troublemaker. I haven't
# 	investigated further
# case "$CI_JOB_NAME" in
#     *"32-bit freebsd"*) : ;;
#     *) freebsd_packages="$freebsd_packages     hwloc" ;;
# esac
freebsd_packages="$freebsd_packages     gmp"
freebsd_packages="$freebsd_packages     gmake"
freebsd_packages="$freebsd_packages     bash"
freebsd_packages="$freebsd_packages     perl5"
freebsd_packages="$freebsd_packages     python3 py39-sqlite3"

while [ $# -gt 0 ] ; do
    case "$1" in
        coverage|clang|gcc|debug|icc) eval "$1=1";;
        *) echo "$1 -> ???" >&2 ; exit 1;;
    esac
    shift
done

# These variables are set in ci/001-environment.sh
if [ "$coverage" ] ; then
    # remove coverage from this round of package selection because we'll
    # install a specific version via pip instead.
    #     debian_packages="$debian_packages     gcovr"
    #     opensuse_packages="$opensuse_packages gcovr"
    #     fedora_packages="$fedora_packages     gcovr"
    #     centos_packages="$centos_packages     gcovr"
    #     alpine_packages="$alpine_packages     gcovr"
    debian_packages="$debian_packages     python3-pip"
    opensuse_packages="$opensuse_packages python3-pip"
    fedora_packages="$fedora_packages     python3-pip"
    centos_packages="$centos_packages     python3-pip"
    alpine_packages="$alpine_packages     py3-pip"
    # vim is needed because we have a bit of ex scripting...
    debian_packages="$debian_packages     lcov vim-nox"
    opensuse_packages="$opensuse_packages lcov vim"
    fedora_packages="$fedora_packages     lcov vim"
    centos_packages="$centos_packages     lcov vim"
    alpine_packages="$alpine_packages     lcov vim"
    if is_freebsd ; then
        echo "coverage -> not on freebsd" >&2
        freebsd_packages="$freebsd_packages   lcov vim"
        # freebsd has no gcovr at the moment, so it's a no-go for now. not
        # sure we expect much benefit in running coverage tests on fbsd as
        # well anyway
        exit 1
    fi
fi

if [ "$gcc32" ] ; then
    if ! is_debian ; then
        echo "multlib -> only debian (fedora,opensuse:IDK ; alpine:no-go)" >&2
        # didn't even check freebsd
        exit 1
    else
        dpkg --add-architecture i386
    fi
    # note that opensuse has gmp-devel-32bit
    debian_packages="$debian_packages     g++-multilib"
    debian_packages="$debian_packages     curl"
    debian_packages="$debian_packages     lzip"
    debian_packages="$debian_packages     libhwloc-dev:i386"
    # we may consider using a system libfmt at some point.
    # debian_packages="$debian_packages     libfmt-dev:i386"
    # fedora_packages="$fedora_packages     g++"
    # centos_packages="$centos_packages     g++"
    # alpine_packages="$alpine_packages     g++"
fi

# The gcc image actually contains a base g++ installation that is in /usr
if [ "$gcc" ] && ! type g++ > /dev/null 2>&1 ; then
    debian_packages="$debian_packages     g++"
    opensuse_packages="$opensuse_packages gcc gcc-c++"
    fedora_packages="$fedora_packages     g++"
    centos_packages="$centos_packages     gcc-c++"
    alpine_packages="$alpine_packages     g++"
    freebsd_packages="$freebsd_packages   gcc"  # this pulls g++ too
fi

if [ "$clang" ] && ! type clang++ > /dev/null 2>&1 ; then
    debian_packages="$debian_packages     clang"
    opensuse_packages="$opensuse_packages clang"
    fedora_packages="$fedora_packages     clang"
    centos_packages="$centos_packages     clang"
    alpine_packages="$alpine_packages     clang"
    freebsd_packages="$freebsd_packages   llvm"
fi

if [ "$checks" ] ; then
    debian_packages="$debian_packages     xsltproc"
    opensuse_packages="$opensuse_packages libxslt-tools"
    fedora_packages="$fedora_packages     libxslt"
    centos_packages="$centos_packages     libxslt"
    alpine_packages="$alpine_packages     libxslt"
    freebsd_packages="$freebsd_packages   libxslt"
fi

if [ "$coverity" ] ; then
    # nothing special at this point. Note that we've tested it on debian
    # only.
    debian_packages="$debian_packages     curl git"
    opensuse_packages="$opensuse_packages curl git"
    fedora_packages="$fedora_packages     curl git"
    centos_packages="$centos_packages     curl git"
    alpine_packages="$alpine_packages     curl git"
    freebsd_packages="$freebsd_packages   curl git"
fi

if [ "$DOCKER_SCRIPT" ] ; then
    debian_packages="$debian_packages sudo git vim-nox gdb"
    opensuse_packages="$opensuse_packages sudo git vim gdb"
    fedora_packages="$fedora_packages sudo git vim gdb"
    centos_packages="$centos_packages sudo git vim gdb"
    alpine_packages="$alpine_packages sudo git vim gdb"
    freebsd_packages="$freebsd_packages sudo git vim gdb"
fi

if is_debian || is_ubuntu ; then
    if [ -x /usr/local/bin/clang ] && ! [ -x /usr/bin/cc ] ; then
        T=$(mktemp -d /tmp/XXXXXXXX)
        F=$($(dirname $0)/phony-packages/clang-usrlocal-0.0/build.sh "$T")
        dpkg -i "$F"
        rm -rf "$T"
    fi
    # intel repos are frequently out of sync, to a point that makes them
    # barely usable. And anyway we don't care: there's no software that
    # we want to pull from these repos anyway.
    find /etc/apt/sources.list.d/ -type f | xargs -r grep -li intel | xargs -r rm
    DEBIAN_FRONTEND=noninteractive apt-get -y update
    DEBIAN_FRONTEND=noninteractive apt-get -y install $debian_packages
elif is_opensuse ; then
    zypper -n install $opensuse_packages
elif is_fedora ; then
    dnf -y install $fedora_packages
elif is_centos ; then
    dnf -y install $centos_packages
elif is_alpine ; then
    # hwloc-dev still in alpine testing.
    cat >> /etc/apk/repositories <<EOF
http://dl-cdn.alpinelinux.org/alpine/edge/testing
EOF
    # is the community repo useful ?
    #http://dl-cdn.alpinelinux.org/alpine/edge/community
    apk update
    apk add $alpine_packages
elif is_freebsd ; then
    env ASSUME_ALWAYS_YES=yes pkg install $freebsd_packages
    if [ "$gcc" ] ; then
        # oh, it's really ugly.
        # with gcc, there's a version clash between the _system_
        # /lib/libgcc_s.so, and the one of the gcc port.
        # https://forums.freebsd.org/threads/freebsd-11-2-libgcc_s-so-1-error.67031/
        # https://wiki.freebsd.org/Ports/libgcc_linking_problem
        mkdir /usr/local/etc/libmap.d
        find /usr/local/lib/gcc* -name '*.so' -o -name '*.so.[0-9]*' | while read xx ; do if [ -e "/lib/$(basename $xx)" ] ; then echo "$(basename "$xx") $xx" ; fi ; done > /usr/local/etc/libmap.d/gcc.conf
    fi
else
    echo "This system is not recognized by our scripts." >&2
    exit 1
fi

if [ "$coverage" ] ; then
    python3 -u -m pip install gcovr==5.0
fi

if [ "$gcc32" ] ; then
    NCPUS=`"$(dirname $0)/utilities/ncpus.sh"`
    cd /tmp/
    curl -O https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz
    tar xf gmp-6.2.1.tar.lz
    cd gmp-6.2.1
    # $GMP is set in ci/001-environment.sh
    ./configure --prefix=$GMP ABI=32
    make -j$NCPUS
    make install
fi

leave_section
