add_subdirectory(googletest-master)
aux_source_directory(${CMAKE_SOURCE_DIR}/tests SRC_LIST)
INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}/googletest-master/googletest/include
	${CMAKE_SOURCE_DIR}/include)
ADD_DEFINITIONS(-std=c++14 -g)
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
add_executable(myTest ${SRC_LIST})
target_link_libraries(myTest myLibrary gtest gtest_main)
add_test(NAME myTest COMMAND myTest)
