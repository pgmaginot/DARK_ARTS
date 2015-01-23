# ########### Full dark_arts testing 
add_test(Linear_Space_Constant_Time_Test ${PROJECT_SOURCE_DIR}/build/dark_arts ${PROJECT_SOURCE_DIR}/testing/TestingInputs/LinearSpaceConstantTime.xml )

set_tests_properties( Linear_Space_Constant_Time_Test PROPERTIES LABELS "dark_arts;MMS")