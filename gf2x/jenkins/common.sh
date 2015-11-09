unset `locale | cut -d= -f1`
export LANG=C

for f in autoconf automake libtool ; do
    if [ -d $HOME/Packages/$m ] ; then
        export PATH=$HOME/Packages/$f/bin:$PATH
    fi
done
# Note: for libtool, we also have to add its aclocal search path to the
# automake aclocal search path (which is
# [automake_prefix]/share/aclocal/dirlist); this is done on all ci nodes.
# 
# For example:
#
#  ls -d /usr/share/aclocal >> ~/Packages/automake-1.15/share/aclocal/dirlist
#
# although, almost by definition, your mileage may vary.


# openbsd has wrappers which strictly require the AUTOCONF_VERSION and
# AUTOMAKE_VERSION variables be set. For most other nodes, this is
# irrelevant and ignored.
export AUTOCONF_VERSION=2.69
if type aclocal-1.15 > /dev/null 2>&1 ; then
    export AUTOMAKE_VERSION=1.15
elif type aclocal-1.14 > /dev/null 2>&1 ; then
    export AUTOMAKE_VERSION=1.14
elif type aclocal-1.13 > /dev/null 2>&1 ; then
    export AUTOMAKE_VERSION=1.13
elif type aclocal > /dev/null ; then
    :
else
    echo "automake (aclocal) not found" >&2
    exit 1
fi

# By default, we do check the FFT interface, but as an option, we also
# allow not checking it.
if ! [ "$DISABLE_FFT" ] ; then
    configure_extra="$configure_extra --enable-fft-interface"
fi

if [ -d "$HOME/Packages/gmp-6.0.0" ] ; then
    configure_extra="$configure_extra --with-gmp=$HOME/Packages/gmp-6.0.0"
elif [ -f "/usr/local/include/gmp.h" ] ; then
    configure_extra="$configure_extra --with-gmp=/usr/local"
fi
