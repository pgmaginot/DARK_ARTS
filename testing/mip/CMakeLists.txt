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