
MAKE_EXTRA_ARGS=

tweak_tree_before_configure() {
    # This is really specific to gf2x
    if [ "$remove_gpl_sources" ] ; then
        ln -sf "toom-gpl-placeholder.c" "$source_tree"/toom-gpl.c
        $ECHO_E "${CSI_BLUE}Removing toom-gpl.c (trying an LGPL test)${CSI_RESET}"
    fi

    if ! [ "$documentation" ] ; then
        MAKE_EXTRA_ARGS="$MAKE_EXTRA_ARGS MAKEINFO=true"
    fi
}

step_configure() {
    if [ "$GMP" ] ; then
        step_configure_autotools_default --with-gmp="$GMP"
    else
        step_configure_autotools_default
    fi
}

build_steps="build"
build_step_name_build="Building"
step_build() {
    (cd "$build_tree" ; "${MAKE}" -j$NCPUS $MAKE_EXTRA_ARGS)
}


check_environment() { : ; }

step_coverage_base() { fatal_error "coverage not set up yet" ; }
coverage_expunge_paths=""
step_coverage_app() { fatal_error "coverage not set up yet" ; }

step_check() {
    # is this "if" really needed after all ?
    if [ "$out_of_source" ] ; then
        (cd "$build_tree" ; make -j$NCPUS check $MAKE_EXTRA_ARGS)
    else
        "${MAKE}" -j$NCPUS check $MAKE_EXTRA_ARGS
    fi
    rc=$?
}

step_coverage_more_artifacts() { : ; }

step_doc() { step_doc_autotools_default "$@" ; }
