set(EXE_NAME Time_Marcher_Initialization)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Time_Marcher_Initialization_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/LinearSpaceConstantTime.xml )

set_tests_properties( Time_Marcher_Initialization_Test PROPERTIES LABELS "MMS;Time_Marcher")
