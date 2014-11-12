
set(EXE_NAME P2_Gauss)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(P2_Gauss_Test ${EXE_NAME}  )

# ###############################################################

set(EXE_NAME P4_NewtonCotes)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(P4_NewtonCotes_Test ${EXE_NAME}  )
set_tests_properties( P4_NewtonCotes_Test PROPERTIES LABELS "quadrature_fem_math; NewtonCotes")

# ###############################################################

set(EXE_NAME P3_Lobatto)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(P3_Lobatto_Test ${EXE_NAME}  )
set_tests_properties( P3_Lobatto_Test PROPERTIES LABELS "quadrature_fem_math; Lobatto")

# ###############################################################

set(EXE_NAME Legendre_Polynomial)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Legendre_Polynomial_Test ${EXE_NAME}  )
set_tests_properties( Legendre_Polynomial_Test PROPERTIES LABELS "quadrature_fem_math; Legendre")

# ###############################################################

set(EXE_NAME Angle_Quadrature)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Angle_Quadrature_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/One_Cell_Angle_Quadrature.xml )
set_tests_properties( Angle_Quadrature_Test PROPERTIES LABELS "quadrature_fem_math; Angular")

# ###############################################################

set(EXE_NAME Crazy_Fem_Strategy)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(Crazy_Fem_Strategy_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/Crazy_Fem_Strategy.xml )
set_tests_properties( Crazy_Fem_Strategy_Test PROPERTIES LABELS "quadrature_fem_math; Fem")

# ###############################################################

# set(EXE_NAME Exact_Integration)

#add_executable(${EXE_NAME} ${EXE_NAME}.cc)
#target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

#add_test(Exact_Integration_Test ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/One_Cell_Angle_Quadrature.xml )
