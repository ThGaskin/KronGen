if (DOXYGEN_FOUND)
    # set the input directories (relative to the top-level source dir)
    set(doxygen_dirs src)

    # configure doxygen
    set(DOXYGEN_QUIET YES)

    set(DOXYGEN_EXTRACT_ALL YES)
    set(DOXYGEN_EXTRACT_PRIVATE YES)
    set(DOXYGEN_EXTRACT_STATIC YES)

    set(DOXYGEN_GENERATE_XML YES)

    set(DOXYGEN_WARNINGS YES)
    set(DOXYGEN_WARN_IF_UNDOCUMENTED YES)
    set(DOXYGEN_WARN_IF_DOC_ERROR YES)
    set(DOXYGEN_WARN_NO_PARAMDOC YES)
    set(DOXYGEN_WARN_LOGFILE ${CMAKE_CURRENT_BINARY_DIR}/doxygen_warnings.log)

    # register the target
    doxygen_add_docs(doxygen_utopia_models
        ${doxygen_dirs}
        # use top-level dir as base directory
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Creating doxygen documentation. Check doxygen/doxygen_warnings.log for emitted warnings."
    )
endif ()
