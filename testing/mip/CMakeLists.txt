set(EXE_NAME MIP_SLXS_Lobatto_Matrices)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(MIP_SLXS_Lobatto_Matrices_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/2Region_Equal_Spacing_MIP_Test.xml )

set_tests_properties( MIP_SLXS_Lobatto_Matrices_Test PROPERTIES LABELS "mip;Lumping;Lobatto")

# ################################################################
