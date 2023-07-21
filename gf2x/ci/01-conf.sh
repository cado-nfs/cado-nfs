#!/usr/bin/env bash

# This is actually __ONLY__ called from ci/40-testsuite.sh now

if ! [ "$sourced_from_testsuite" ] ; then
    . "${CI_PATH:-$(dirname $0)}/000-functions.sh"
    . "${CI_PATH:-$(dirname $0)}/001-environment.sh"
fi

enter_section configuration Configuring
autoreconf -i

if [ "$out_of_source" ] ; then
    if ! [ "$build_tree" ] ; then
        build_tree="${TMPDIR:-/tmp}/$CI_JOB_NAME"
        # spaces in dir names don't work, mostly because of libtool
        # (look at gf2x/fft/libgf2x-fft.la)
        # This substitution is bash-only, but this should be fine to 
        # have in a conditional that non-bash skips over
        build_tree="${build_tree// /_}"
    fi
    if ! [ -d "$build_tree" ] ; then
        mkdir "$build_tree"
    fi
    $ECHO_E "${CSI_BLUE}out-of-source build tree is $build_tree${CSI_RESET}"
    configure_call="$source_tree/configure"
else
    build_tree=$PWD
    configure_call="./configure"
fi

export build_tree
export source_tree
export configure_call

if [ "$remove_gpl_sources" ] ; then
    ln -sf "toom-gpl-placeholder.c" "$source_tree"/toom-gpl.c
    $ECHO_E "${CSI_BLUE}Removing toom-gpl.c (trying an LGPL test)${CSI_RESET}"
fi

(cd "$build_tree" ; "$configure_call")

if [ "$do_make_dist" ] ; then
    eval $(cd "$build_tree" > /dev/null ; grep '^\(PACKAGE_\(TARNAME\|VERSION\)\)' Makefile | tr -d \  )
    $ECHO_E "${CSI_BLUE}Creating archive ${PACKAGE_TARNAME}-${PACKAGE_VERSION}.tar.gz${CSI_RESET}"
    (cd "$build_tree" ; "${MAKE}" dist)
    cd "${TMPDIR:-/tmp}"
    source_tree="$PWD/${PACKAGE_TARNAME}-${PACKAGE_VERSION}"
    $ECHO_E "${CSI_BLUE}Unpacking in $source_tree, and configuring there${CSI_RESET}"
    tar xf "$build_tree/${PACKAGE_TARNAME}-${PACKAGE_VERSION}.tar.gz"
    # hack, and change source tree and build tree.
    if [ "$out_of_source" ] ; then
        rm -rf "$build_tree"
        build_tree="${TMPDIR:-/tmp}/$CI_JOB_NAME"
        build_tree="${build_tree// /_}"
        mkdir "$build_tree"
        $ECHO_E "${CSI_BLUE}out-of-source build tree is $build_tree${CSI_RESET}"
        configure_call="$source_tree/configure"
    else
        $ECHO_E "${CSI_BLUE}Changing source directory to $source_tree${CSI_RESET}"
        cd $source_tree
        build_tree=$PWD
    fi
    export build_tree
    export source_tree
    (cd "$build_tree" ; "$configure_call")
fi

leave_section

