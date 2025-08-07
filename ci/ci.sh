
export project_name=cado-nfs
export build_system=cmake

needs_bc=1
needs_python=1
needs_perl=1
needs_optional_hwloc=1
needs_optional_ecm=
needs_optional_fmt=1
needs_gmp=1

if is_debian || is_ubuntu || is_fedora ; then
    # it's easy when we have a package.
    needs_optional_ecm=1
    # see also special case for "merge coverage tests" further down
fi

case "$JOB_NAME" in
    *"under valgrind"*)
        valgrind=1
        export VALGRIND=1
        export USE_ONLY_ASSEMBLY_INSTRUCTIONS_THAT_VALGRIND_KNOWS_ABOUT=1
        export TIMEOUT_SCALE=10
        case "$JOB_NAME" in
            *32-bit*)
                echo "valgrind testing is practically hopeless under 32-bit"
                echo "See https://bugs.kde.org/show_bug.cgi?id=337475"
                echo
                echo "We're making this test artificially succeed"
                exit 0
                # Things like shlx (in the BMI2 set) are valid in 32-bit
                # mode, and can definitely be emitted by a compiler, or
                # can be present in libraries such as gmp. We have no way
                # to guarantee that no such instructions are encountered,
                # and valgrind is unwilling to ramp up support for these
                # instructions (which is understandable)
                ;;
        esac
        ;;
    *"coverage tests"*)
        # With some tests, the coverage test time goes to the roof, but
        # it's not always so. A blanket TIMEOUT_SCALE is probably a gross
        # fix, but we can live with it.
        export TIMEOUT_SCALE=2

        # # It's probably debatable. If we _do_ get more coverage with
        # # expensive checks, then frankly, I would consider it a bit of a
        # # bug: it would be better to get the same coverage with the
        # # normal tests.
        # export CHECKS_EXPENSIVE=1

        debian_packages="$debian_packages     lcov libmpc-dev libmpfr-dev"
        alpine_packages="$alpine_packages     lcov mpc1-dev mpfr-dev"
        needs_optional_ecm=1
        # see also special case for "merge coverage tests" further down
        if ! is_debian && ! is_alpine ; then
            echo "lcov: only on debian|alpine" >&2
            # because I'm lazy, and also I'm not sure there would be a point
            # in doing it on several systems anyway.
            exit 1
        fi
        ;;
    *"using package libfmt-dev"*)
        export install_package_libfmt_dev=1
        ;;
    *"mysql specific"*)
        export mysql=1
        ;;
esac

case "$JOB_NAME" in
    *"merge coverage tests"*)
        # here we're not even tied to a compiler in particular, so let's
        # avoid requesting software that requires recompilation.
        needs_optional_ecm=
    ;;
esac

project_package_selection() {
    if [ "$valgrind" ] ; then
        echo " + valgrind is set"
        debian_packages="$debian_packages     valgrind"
        if ! is_debian ; then
            echo "valgrind: only on debian" >&2
            # because I'm lazy, and also I'm not sure there would be a point
            # in doing it on several systems anyway.
            exit 1
        fi
    fi

    if [ "$mysql" ] ; then
        echo " + mysql is set"
        # https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=923347
        pip_packages="$pip_packages mysql.connector"
        debian_packages="$debian_packages     mariadb-server"
    fi

    if [ "$install_package_libfmt_dev" ] ; then
        echo " + libfmt_dev is set"
        debian_packages="$debian_packages     libfmt-dev"
        if ! is_debian && ! is_ubuntu ; then
            echo "libfmt-dev: only on debian" >&2
            # because I'm lazy, and also I'm not sure there would be a point
            # in doing it on several systems anyway.
            exit 1
        fi
    fi
   
    # bzip2 is used in the tests, but after all it makes sense that we
    # pay attention to our good support of some compression tool not
    # being available...
    #
    # debian_packages="$debian_packages     bzip2"
    # opensuse_packages="$opensuse_packages bzip2"
    # fedora_packages="$fedora_packages     bzip2"
    # centos_packages="$centos_packages     bzip2"
    # alpine_packages="$alpine_packages     bzip2"

    if [ "$needs_python" ] ; then
        echo " + needs_python is set"
        debian_packages="$debian_packages     python3-flask python3-requests"
        opensuse_packages="$opensuse_packages python3-Flask python3-requests"
        fedora_packages="$fedora_packages     python3-flask python3-requests openssl"
        centos_packages="$centos_packages     python3-requests"
        if is_centos ; then
            # centos has no python3-flask package, but we can install it
            # via pip
            pip_packages="$pip_packages flask"
        fi
        alpine_packages="$alpine_packages     py3-flask py3-requests"
        # py311-sqlite3 is in the python stdlib, but trimmed on on fbsd
        freebsd_packages="$freebsd_packages   py311-sqlite3 py311-flask py311-requests"
        freebsd_packages="$freebsd_packages   pkgconf" # added PZ
    fi

    if [ "$needs_optional_ecm" ] ; then
        # ecm is heasy when we have a distro package, which is the case
        # with debian/ubuntu/fedora. We also want to make sure that we
        # have it even with alpine linux, which requires some
        # recompiling.  (see the function after_package_install() down
        # below). Building gmp-ecm requires m4.
        echo " + needs_optional_ecm is set"
        debian_packages="$debian_packages     libecm-dev"
        fedora_packages="$fedora_packages     gmp-ecm-devel"
        alpine_packages="$alpine_packages     curl m4"
    fi

    # add this so that we get the gdb tests as well (at least with the
    # shared libs on debian-testing case)
    debian_packages="$debian_packages     gdb"
}

after_package_install() {
    # This is still run as root.
    if [ "$mysql" ] ; then
        /etc/init.d/mariadb start
        # ci jobs will run all this as root, so we don't really care. But
        # for runs with ci/ci/debug.sh, we do.
        mysql -e "CREATE USER IF NOT EXISTS 'hostuser'@localhost IDENTIFIED VIA unix_socket;"
        mysql -e "GRANT ALL PRIVILEGES ON cado_nfs.* TO hostuser@localhost IDENTIFIED VIA unix_socket;"
    fi
    mkdir -p /etc/gdb
    echo "set auto-load safe-path /" > /etc/gdb/gdbinit

    if [ "$needs_optional_ecm" ] && ! (is_debian || is_ubuntu || is_fedora) ; then
        url=https://gitlab.inria.fr/-/project/24244/uploads/ad3e5019fef98819ceae58b78f4cce93/ecm-7.0.6.tar.gz
        (
            NCPUS=`"${CI_PATH:-$(dirname $0)}/utilities/ncpus.sh"`
            # we don't _really_ need to cd. But OTOH coverage tends to
            # aggressively index code under $PWD, so it's better to keep
            # this one out.
            cd /tmp/
            # $needs_optional_ecm pulls m4 under alpine. As for other
            # systems, I don't know... For now this piece of code is not
            # triggered on other systems, I think.
            # apk update
            # apk add gcc m4 gmp-dev musl-dev make
            curl -O "$url"
            tar xf $(basename "$url")
            cd ecm-7.0.6
            ./configure
            make -j$NCPUS
            make install
        )
    fi
}

# Note: most of the interesting stuff is in ci.bash
