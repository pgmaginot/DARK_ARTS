
set(EXE_NAME P2_Gauss)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(P2_Gauss_Test ${EXE_NAME}  )

# ###############################################################

set(EXE_NAME P4_NewtonCotes)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(P4_NewtonCotes_Test ${EXE_NAME}  )

# ###############################################################

set(EXE_NAME P3_Lobatto)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(P3_Lobatto_Test ${EXE_NAME}  )
