
set(EXE_NAME LogSpacing)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(LogSpacingTest ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/2Region_Log_Spacing.xml )

#set(EXE_NAME testRosanov)

#add_executable(${EXE_NAME} ${EXE_NAME}.cc)
#target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

#add_test(TestRosanov ${EXE_NAME} )