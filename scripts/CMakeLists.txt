
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
    api_server.py
    cadofactor_tools/__init__.py
    cadofactor_tools/certificate.py
    cadofactor_tools/upload_dir.py
    database/__init__.py
    database/access.py
    database/base.py
    database/dictdb.py
    database/factory.py
    database/misc.py
    database/mysql.py
    database/sqlite3.py
    database/table.py
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
        get_filename_component(d ${f} DIRECTORY)
        install(FILES ${m}/${f} DESTINATION ${LIBSUFFIX}/scripts/${m}/${d})
    endforeach()
endforeach()

configure_file(descent.py descent.py COPYONLY)
foreach(m descent)
    foreach(f ${descent_scripts})
        configure_file(${m}/${f} ${m}/${f} COPYONLY)
    endforeach()
endforeach()

install(PROGRAMS descent.py DESTINATION ${LIBSUFFIX}/scripts)
