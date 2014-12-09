set(EXE_NAME SLXS_MMS_Sweep_Matrix)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(SLXS_MMS_Sweep_Matrix_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/SLXS_MMS_Sweep_Matrix.xml )

set_tests_properties( SLXS_MMS_Sweep_Matrix_Test PROPERTIES LABELS "MMS;transport_sweep;Lumping;Lobatto")

