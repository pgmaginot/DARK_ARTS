
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

# #######################################################################################

set(EXE_NAME Interpolating)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Interpolating_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/Crazy_Fem_Strategy.xml )

set_tests_properties( Interpolating_Test PROPERTIES LABELS "materials")

# #######################################################################################

set(EXE_NAME SLXS)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(SLXS_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/2Region_Log_Spacing.xml )

set_tests_properties( SLXS_Test PROPERTIES LABELS "materials;SLXS")

# #######################################################################################

set(EXE_NAME MMS_Sources)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MMS_Sources_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/MMS_Source_Test.xml )

set_tests_properties( MMS_Sources_Test PROPERTIES LABELS "materials;MMS")

# #######################################################################################

set(EXE_NAME MMS_Components)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MMS_Components ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/MMS_Source_Test.xml )

set_tests_properties( MMS_Components PROPERTIES LABELS "materials;MMS")