#
#  Need to test all of the following:  V_Intensity_Update, Intensity_Update_Grey , V_Temperature_Update, Temperature_Update_Grey , WGRS_FP_Sweeps, V_WGRS
#
#  
#

set(EXE_NAME MMS_Intensity_Update)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MMS_Intensity_Update_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/LinearSpaceConstantTime.xml )

set_tests_properties( MMS_Intensity_Update_Test PROPERTIES LABELS "MMS;Intensity_Update")

########################################################

set(EXE_NAME FP_Sweeps)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(FP_Sweeps_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/LinearSpaceConstantTime.xml )

set_tests_properties( FP_Sweeps_Test PROPERTIES LABELS "MMS;WGRS;transport_sweep")

########################################################

set(EXE_NAME Single_Thermal_Iteration)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Single_Thermal_Iteration_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/LinearSpaceConstantTime.xml )

set_tests_properties( Single_Thermal_Iteration_Test PROPERTIES LABELS "MMS;WGRS;transport_sweep;Intensity_Update;Temperature_Update")