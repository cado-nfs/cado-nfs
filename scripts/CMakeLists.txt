
set(cadofactor_scripts
    cadocommand.py
    cadologger.py
    cadoparams.py
    cadoprograms.py
    cadotask.py
    patterns.py
    upload.py
    workunit.py
    toplevel.py
    wuserver.py
    # signalhandler.py
    # verify.py
    cadofactor_tools/__init__.py
    cadofactor_tools/certificate.py
)

set(descent_scripts
    __main__.py
    __init__.py
    descent_helper_asyncio.py
    descent_helper_fallback.py
    descent_lower_class.py
    descent_middle_class.py
    descent_upper_class.py
    descent_utils.py
    general_class.py
    ideals_above_p.py
    log_base.py
)

foreach(m cadofactor descent)
    foreach(f ${${m}_scripts})
        install(FILES ${m}/${f} DESTINATION ${LIBSUFFIX}/scripts/${m})
    endforeach()
endforeach()

# This one is a bit special since we do @ replacement there.
configure_file(cadofactor/wudb.py cadofactor/wudb.py ESCAPE_QUOTES @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/cadofactor/wudb.py DESTINATION ${LIBSUFFIX}/scripts/cadofactor)

configure_file(descent.py descent.py COPYONLY)
foreach(m descent)
    foreach(f ${descent_scripts})
        configure_file(${m}/${f} ${m}/${f} COPYONLY)
    endforeach()
endforeach()

install(PROGRAMS descent.py DESTINATION ${LIBSUFFIX}/scripts)
