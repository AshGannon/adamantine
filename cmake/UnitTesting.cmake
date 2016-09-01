function(adamantine_ADD_BOOST_TEST TEST_NAME)
    add_executable(${TEST_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.cc ${tests_SOURCES})
    target_include_directories(${TEST_NAME} SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
    target_link_libraries(${TEST_NAME} PUBLIC ${Boost_CHRONO_LIBRARY})
    target_link_libraries(${TEST_NAME} PUBLIC ${Boost_FILESYSTEM_LIBRARY})
    target_link_libraries(${TEST_NAME} PUBLIC ${Boost_MPI_LIBRARY})
    target_link_libraries(${TEST_NAME} PUBLIC ${Boost_PROGRAM_OPTIONS_LIBRARY})
    target_link_libraries(${TEST_NAME} PUBLIC ${Boost_TIMER_LIBRARY})
    target_link_libraries(${TEST_NAME} PUBLIC ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    target_include_directories(${TEST_NAME} SYSTEM PUBLIC ${DEAL_II_INCLUDE_DIRS})
    target_link_libraries(${TEST_NAME} PUBLIC ${DEAL_II_LIBRARIES})
    target_link_libraries(${TEST_NAME} LINK_PUBLIC Adamantine)
    set_target_properties(${TEST_NAME} PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED ON
    )
    if(ARGN)
        set(NUMBER_OF_PROCESSES_TO_EXECUTE ${ARGN})
    else()
        set(NUMBER_OF_PROCESSES_TO_EXECUTE 1)
    endif()
    foreach(NPROC ${NUMBER_OF_PROCESSES_TO_EXECUTE})
        add_test(
            NAME ${TEST_NAME}_cpp_${NPROC}
            COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} ${CMAKE_BINARY_DIR}/bin/${TEST_NAME}
        )
        set_tests_properties(${TEST_NAME}_cpp_${NPROC} PROPERTIES
            PROCESSORS ${NPROC}
        )
    endforeach()
endfunction()