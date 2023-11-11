
export project_name=cado-nfs
: ${CI_PROJECT_URL=https://gitlab.inria.fr/cado-nfs/cado-nfs}
export build_system=cmake

needs_bc=1
needs_python=1
needs_perl=1
needs_optional_hwloc=1
needs_optional_ecm=1
needs_optional_fmt=1
needs_gmp=1

case "$CI_JOB_NAME" in
    *"under valgrind"*)
        valgrind=1
        export USE_ONLY_ASSEMBLY_INSTRUCTIONS_THAT_VALGRIND_KNOWS_ABOUT=1
        case "$CI_JOB_NAME" in
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
        # opensuse_packages="$opensuse_packages python3-pip"
        # fedora_packages="$fedora_packages     python3-pip"
        # centos_packages="$centos_packages     python3-pip"
        # alpine_packages="$alpine_packages     py3-pip"
    fi
}

# Note: most of the interesting stuff is in ci.bash
