# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ian/Engineer/colors

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ian/Engineer/colors/build

# Include any dependencies generated for this target.
include source/CMakeFiles/SRC_LIST.dir/depend.make

# Include the progress variables for this target.
include source/CMakeFiles/SRC_LIST.dir/progress.make

# Include the compile flags for this target's objects.
include source/CMakeFiles/SRC_LIST.dir/flags.make

source/CMakeFiles/SRC_LIST.dir/requires:
.PHONY : source/CMakeFiles/SRC_LIST.dir/requires

source/CMakeFiles/SRC_LIST.dir/clean:
	cd /home/ian/Engineer/colors/build/source && $(CMAKE_COMMAND) -P CMakeFiles/SRC_LIST.dir/cmake_clean.cmake
.PHONY : source/CMakeFiles/SRC_LIST.dir/clean

source/CMakeFiles/SRC_LIST.dir/depend:
	cd /home/ian/Engineer/colors/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ian/Engineer/colors /home/ian/Engineer/colors/source /home/ian/Engineer/colors/build /home/ian/Engineer/colors/build/source /home/ian/Engineer/colors/build/source/CMakeFiles/SRC_LIST.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : source/CMakeFiles/SRC_LIST.dir/depend

