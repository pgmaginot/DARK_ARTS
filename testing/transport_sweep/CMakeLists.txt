set(EXE_NAME SLXS_MMS_Sweep_Matrix)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(SLXS_MMS_Sweep_Matrix_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/SLXS_MMS_Sweep_Matrix.xml )

set_tests_properties( SLXS_MMS_Sweep_Matrix_Test PROPERTIES LABELS "MMS;transport_sweep;Lumping;Lobatto")

# ################################################################

# set(EXE_NAME Simple_BC)

# add_executable(${EXE_NAME} ${EXE_NAME}.cc)
# target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

# add_test(Simple_BC_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/SLXS_MMS_Sweep_Matrix.xml )

# set_tests_properties( Simple_BC_Test PROPERTIES LABELS "transport_sweep;BC")

# ################################################################


set(EXE_NAME MMS_BC)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MMS_BC_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/MMS_Source_Test.xml )

set_tests_properties( MMS_BC_Test PROPERTIES LABELS "MMS;transport_sweep;BC")

# ################################################################

set(EXE_NAME ZeroChangeSweep)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(ZeroChangeSweep_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/LinearSpaceConstantTime.xml )

set_tests_properties( ZeroChangeSweep_Test PROPERTIES LABELS "transport_sweep")

# ################################################################

set(EXE_NAME Upwinding)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Upwinding_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/MMS_Source_Test.xml )

set_tests_properties( Upwinding_Test PROPERTIES LABELS "MMS;transport_sweep;Upwind")