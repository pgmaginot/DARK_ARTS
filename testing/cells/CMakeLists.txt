
set(EXE_NAME LogSpacing)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(LogSpacingTest ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/2Region_Log_Spacing.xml )

set_tests_properties( LogSpacingTest PROPERTIES LABELS "cells;Log")

# ###################################################################

set(EXE_NAME EqualSpacing)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(EqualSpacingTest ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/1Region_Equal_Spacing.xml )

set_tests_properties( EqualSpacingTest PROPERTIES LABELS "cells;Equal")

# ###################################################################

set(EXE_NAME LogSpacingRefinement)

add_executable(${EXE_NAME} ${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${${PROJECT_NAME}_LIB_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})

add_test(LogSpacingRefinementTest ${EXE_NAME} ${PROJECT_SOURCE_DIR}/testing/TestingInputs/2Region_Log_Spacing_Refinement.xml )

set_tests_properties( LogSpacingRefinementTest PROPERTIES LABELS "cells;Log;refinement")