
set(EXE_NAME Unity_Grey_Planck)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Unity_Grey_Planck_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/1Region_Equal_Spacing.xml )

set_tests_properties( Unity_Grey_Planck_Test PROPERTIES LABELS "materials;Grey_Planck")

# #######################################################################################

set(EXE_NAME Materials_Interface)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Materials_Interface_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/One_Cell_Angle_Quadrature.xml )

set_tests_properties( Materials_Interface_Test PROPERTIES LABELS "materials")