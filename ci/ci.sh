
export project_name=cado-nfs
export build_system=cmake

needs_bc=1
needs_python=1
needs_perl=1
needs_optional_hwloc=1
needs_optional_ecm=1
needs_optional_fmt=1
needs_gmp=1

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

        debian_packages="$debian_packages     lcov"
        alpine_packages="$alpine_packages     lcov"
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
}

# Note: most of the interesting stuff is in ci.bash
