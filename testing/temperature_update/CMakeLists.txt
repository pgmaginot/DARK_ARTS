
set(EXE_NAME Grey_Update_Matrices)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Grey_Update_Matrices_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/LinearSpaceLinearTime.xml )

set_tests_properties( Grey_Update_Matrices_Test PROPERTIES LABELS "MMS;Temperature_Update;Temperature_Matrix")