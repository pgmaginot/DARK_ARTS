
set(EXE_NAME Ramp_Time_Data)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Ramp_Time_Data_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/2Region_Log_Spacing.xml )

set_tests_properties( Ramp_Time_Data_Test PROPERTIES LABELS "time_data;Ramp")

#####################################################################

set(EXE_NAME Vector_Time_Data)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Vector_Time_Data_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/1Region_Equal_Spacing.xml )

set_tests_properties( Vector_Time_Data_Test PROPERTIES LABELS "time_data;Vector")

# #####################################################################

set(EXE_NAME Exponential_Time_Data)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Exponential_Time_Data_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/Exponential_Time.xml )

set_tests_properties( Exponential_Time_Data_Test PROPERTIES LABELS "time_data;Exponential")


#####################################################################

set(EXE_NAME SDIRK_22)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(SDIRK_22_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/SDIRK_22.xml )

set_tests_properties( SDIRK_22_Test PROPERTIES LABELS "time_data")


#####################################################################
set(EXE_NAME SDIRK_33)


add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(SDIRK_33_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/SDIRK_33.xml )

set_tests_properties( SDIRK_33_Test PROPERTIES LABELS "time_data")

#####################################################################
set(EXE_NAME MMS_Test)


add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MMS_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/SDIRK_22.xml )

set_tests_properties( MMS_Test PROPERTIES LABELS "time_data")