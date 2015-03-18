set(EXE_NAME MIP_SLXS_Lobatto_Matrices)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MIP_SLXS_Lobatto_Matrices_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/2Region_Equal_Spacing_MIP_Test.xml )

set_tests_properties( MIP_SLXS_Lobatto_Matrices_Test PROPERTIES LABELS "mip;Lumping;Lobatto")

# ################################################################

set(EXE_NAME Kappa_Calculator)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Kappa_Calculator_Test ${EXE_NAME})

set_tests_properties( Kappa_Calculator_Test PROPERTIES LABELS "mip")

# ################################################################

set(EXE_NAME MIP_Local_Matrix_Assembler)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MIP_Local_Matrix_Assembler_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/MIP_Assembler_Test.xml )

set_tests_properties( MIP_Local_Matrix_Assembler_Test PROPERTIES LABELS "mip")

# ################################################################

set(EXE_NAME MIP_Exact_Matrices)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MIP_Exact_Matrices_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/2Region_Equal_Spacing_MIP_Test.xml )

set_tests_properties( MIP_Exact_Matrices_Test PROPERTIES LABELS "mip;Lumping;Lobatto")

# ################################################################

set(EXE_NAME Row_to_Col_Major_Results)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Row_to_Col_Major_Results_Test ${EXE_NAME}  )

set_tests_properties( Row_to_Col_Major_Results_Test PROPERTIES LABELS "mip;Eigen")

# ################################################################


set(EXE_NAME MIP_Arbitrary_R_Mat_Check)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MIP_Arbitrary_R_Mat_Check ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/MMS1_TL_P1_4_Cells_DSA.xml  )

set_tests_properties( MIP_Arbitrary_R_Mat_Check PROPERTIES LABELS "mip")

# ################################################################


set(EXE_NAME MIP_Matrix_Assembly)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MIP_Matrix_Assembly ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/Constant_In_Space_MIP_Check.xml  )

set_tests_properties( MIP_Matrix_Assembly PROPERTIES LABELS "mip")

# ################################################################