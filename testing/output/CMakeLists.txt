set(EXE_NAME No_Difference_L2)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(No_Difference_L2_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/LinearSpaceConstantTime.xml )

set_tests_properties( No_Difference_L2_Test PROPERTIES LABELS "MMS;Error;Output")

# ##############################################################

set(EXE_NAME Linear_Difference_L2)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Linear_Difference_L2_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/QuadraticSpaceConstantTime.xml )

set_tests_properties( Linear_Difference_L2_Test PROPERTIES LABELS "MMS;Error;Output")