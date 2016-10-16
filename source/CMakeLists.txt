aux_source_directory(${CMAKE_SOURCE_DIR}/source SRC_LIST)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
ADD_DEFINITIONS(-std=c++14 -g)
list (REMOVE_ITEM SRC_LIST main.cpp ini)
add_library (myLibrary ${SRC_LIST})
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
add_executable(${PROJECT_NAME}.exe main.cpp)
target_link_libraries(${PROJECT_NAME}.exe myLibrary)

