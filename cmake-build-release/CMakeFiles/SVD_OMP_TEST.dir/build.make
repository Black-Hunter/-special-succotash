# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /home/jido/clion-2019.3.3/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/jido/clion-2019.3.3/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/cmake-build-release

# Include any dependencies generated for this target.
include CMakeFiles/SVD_OMP_TEST.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/SVD_OMP_TEST.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SVD_OMP_TEST.dir/flags.make

CMakeFiles/SVD_OMP_TEST.dir/Test/test.cpp.o: CMakeFiles/SVD_OMP_TEST.dir/flags.make
CMakeFiles/SVD_OMP_TEST.dir/Test/test.cpp.o: ../Test/test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/SVD_OMP_TEST.dir/Test/test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SVD_OMP_TEST.dir/Test/test.cpp.o -c /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/Test/test.cpp

CMakeFiles/SVD_OMP_TEST.dir/Test/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SVD_OMP_TEST.dir/Test/test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/Test/test.cpp > CMakeFiles/SVD_OMP_TEST.dir/Test/test.cpp.i

CMakeFiles/SVD_OMP_TEST.dir/Test/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SVD_OMP_TEST.dir/Test/test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/Test/test.cpp -o CMakeFiles/SVD_OMP_TEST.dir/Test/test.cpp.s

CMakeFiles/SVD_OMP_TEST.dir/Main/svd.cpp.o: CMakeFiles/SVD_OMP_TEST.dir/flags.make
CMakeFiles/SVD_OMP_TEST.dir/Main/svd.cpp.o: ../Main/svd.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/SVD_OMP_TEST.dir/Main/svd.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SVD_OMP_TEST.dir/Main/svd.cpp.o -c /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/Main/svd.cpp

CMakeFiles/SVD_OMP_TEST.dir/Main/svd.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SVD_OMP_TEST.dir/Main/svd.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/Main/svd.cpp > CMakeFiles/SVD_OMP_TEST.dir/Main/svd.cpp.i

CMakeFiles/SVD_OMP_TEST.dir/Main/svd.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SVD_OMP_TEST.dir/Main/svd.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/Main/svd.cpp -o CMakeFiles/SVD_OMP_TEST.dir/Main/svd.cpp.s

CMakeFiles/SVD_OMP_TEST.dir/Test/catch_main.cpp.o: CMakeFiles/SVD_OMP_TEST.dir/flags.make
CMakeFiles/SVD_OMP_TEST.dir/Test/catch_main.cpp.o: ../Test/catch_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/SVD_OMP_TEST.dir/Test/catch_main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SVD_OMP_TEST.dir/Test/catch_main.cpp.o -c /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/Test/catch_main.cpp

CMakeFiles/SVD_OMP_TEST.dir/Test/catch_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SVD_OMP_TEST.dir/Test/catch_main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/Test/catch_main.cpp > CMakeFiles/SVD_OMP_TEST.dir/Test/catch_main.cpp.i

CMakeFiles/SVD_OMP_TEST.dir/Test/catch_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SVD_OMP_TEST.dir/Test/catch_main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/Test/catch_main.cpp -o CMakeFiles/SVD_OMP_TEST.dir/Test/catch_main.cpp.s

# Object files for target SVD_OMP_TEST
SVD_OMP_TEST_OBJECTS = \
"CMakeFiles/SVD_OMP_TEST.dir/Test/test.cpp.o" \
"CMakeFiles/SVD_OMP_TEST.dir/Main/svd.cpp.o" \
"CMakeFiles/SVD_OMP_TEST.dir/Test/catch_main.cpp.o"

# External object files for target SVD_OMP_TEST
SVD_OMP_TEST_EXTERNAL_OBJECTS =

SVD_OMP_TEST: CMakeFiles/SVD_OMP_TEST.dir/Test/test.cpp.o
SVD_OMP_TEST: CMakeFiles/SVD_OMP_TEST.dir/Main/svd.cpp.o
SVD_OMP_TEST: CMakeFiles/SVD_OMP_TEST.dir/Test/catch_main.cpp.o
SVD_OMP_TEST: CMakeFiles/SVD_OMP_TEST.dir/build.make
SVD_OMP_TEST: CMakeFiles/SVD_OMP_TEST.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable SVD_OMP_TEST"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SVD_OMP_TEST.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SVD_OMP_TEST.dir/build: SVD_OMP_TEST

.PHONY : CMakeFiles/SVD_OMP_TEST.dir/build

CMakeFiles/SVD_OMP_TEST.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SVD_OMP_TEST.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SVD_OMP_TEST.dir/clean

CMakeFiles/SVD_OMP_TEST.dir/depend:
	cd /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/cmake-build-release /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/cmake-build-release /home/jido/Desktop/Uni_Jena_Informatik_Master/Dritte_Semester/Algorithm_Engnieering/svd_GIT/cmake-build-release/CMakeFiles/SVD_OMP_TEST.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SVD_OMP_TEST.dir/depend

