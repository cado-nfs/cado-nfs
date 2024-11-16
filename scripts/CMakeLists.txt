
set(cadofactor_scripts
    cadocommand.py
    cadologger.py
    cadoparams.py
    cadoprograms.py
    cadotask.py
    patterns.py
    workunit.py
    toplevel.py
    wudb.py
    # signalhandler.py
    # verify.py
    api_server/__init__.py
    api_server/api.py
    api_server/launcher.py
    cadofactor_tools/__init__.py
    cadofactor_tools/certificate.py
    cadofactor_tools/upload_dir.py
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

configure_file(descent.py descent.py COPYONLY)
foreach(m descent)
    foreach(f ${descent_scripts})
        configure_file(${m}/${f} ${m}/${f} COPYONLY)
    endforeach()
endforeach()

install(PROGRAMS descent.py DESTINATION ${LIBSUFFIX}/scripts)
