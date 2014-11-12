
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