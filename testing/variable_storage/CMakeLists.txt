
set(EXE_NAME ConstantTemperature)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(ConstantTemperatureTest ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/2Region_Log_Spacing.xml )

set_tests_properties( ConstantTemperatureTest PROPERTIES LABELS "variable_storage;Temperature")

# ###################################################################

set(EXE_NAME ConstantIntensityIC)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(ConstantIntensityICTest ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/2Region_Log_Spacing.xml )

set_tests_properties( ConstantIntensityICTest PROPERTIES LABELS "variable_storage;Radiation")

# ###################################################################

set(EXE_NAME IntensityMomentData)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(IntensityMomentDataTest ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/2Region_Log_Spacing.xml )

set_tests_properties( IntensityMomentDataTest PROPERTIES LABELS "variable_storage;Radiation")

# ###################################################################

set(EXE_NAME OutputTrial)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(OutputTrialTest ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/1Region_Equal_Spacing.xml )

set_tests_properties( OutputTrialTest PROPERTIES LABELS "output")

# ###################################################################

set(EXE_NAME MMS_Temperature_IC)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MMS_Temperature_IC_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/MMS_Source_Test.xml )

set_tests_properties( MMS_Temperature_IC_Test PROPERTIES LABELS "variable_storage;Temperature;MMS")

# ###################################################################

set(EXE_NAME MMS_Intensity_IC)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MMS_Intensity_IC_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/MMS_Source_Test.xml )

set_tests_properties( MMS_Intensity_IC_Test PROPERTIES LABELS "variable_storage;Radiation;MMS")

# ###################################################################

set(EXE_NAME ErrPhi)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(ErrPhi_Test ${EXE_NAME} )

set_tests_properties( ErrPhi_Test PROPERTIES LABELS "variable_storage;Radiation")

# ###################################################################

set(EXE_NAME ErrTemperature)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(ErrTemperature_Test ${EXE_NAME} )

set_tests_properties( ErrTemperature_Test PROPERTIES LABELS "variable_storage;Temperature")
