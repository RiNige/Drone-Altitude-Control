# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.16.4/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.16.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/tianwang/dropbox/HPC/coursework/Code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/tianwang/dropbox/HPC/coursework/Code/build

# Include any dependencies generated for this target.
include CMakeFiles/myprog.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/myprog.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/myprog.dir/flags.make

CMakeFiles/myprog.dir/LidDrivenCavitySolver.cpp.o: CMakeFiles/myprog.dir/flags.make
CMakeFiles/myprog.dir/LidDrivenCavitySolver.cpp.o: ../LidDrivenCavitySolver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/tianwang/dropbox/HPC/coursework/Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/myprog.dir/LidDrivenCavitySolver.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/myprog.dir/LidDrivenCavitySolver.cpp.o -c /Users/tianwang/dropbox/HPC/coursework/Code/LidDrivenCavitySolver.cpp

CMakeFiles/myprog.dir/LidDrivenCavitySolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/myprog.dir/LidDrivenCavitySolver.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/tianwang/dropbox/HPC/coursework/Code/LidDrivenCavitySolver.cpp > CMakeFiles/myprog.dir/LidDrivenCavitySolver.cpp.i

CMakeFiles/myprog.dir/LidDrivenCavitySolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/myprog.dir/LidDrivenCavitySolver.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/tianwang/dropbox/HPC/coursework/Code/LidDrivenCavitySolver.cpp -o CMakeFiles/myprog.dir/LidDrivenCavitySolver.cpp.s

# Object files for target myprog
myprog_OBJECTS = \
"CMakeFiles/myprog.dir/LidDrivenCavitySolver.cpp.o"

# External object files for target myprog
myprog_EXTERNAL_OBJECTS =

myprog: CMakeFiles/myprog.dir/LidDrivenCavitySolver.cpp.o
myprog: CMakeFiles/myprog.dir/build.make
myprog: /usr/local/Cellar/open-mpi/4.0.2/lib/libmpi.dylib
myprog: CMakeFiles/myprog.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/tianwang/dropbox/HPC/coursework/Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable myprog"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/myprog.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/myprog.dir/build: myprog

.PHONY : CMakeFiles/myprog.dir/build

CMakeFiles/myprog.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/myprog.dir/cmake_clean.cmake
.PHONY : CMakeFiles/myprog.dir/clean

CMakeFiles/myprog.dir/depend:
	cd /Users/tianwang/dropbox/HPC/coursework/Code/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/tianwang/dropbox/HPC/coursework/Code /Users/tianwang/dropbox/HPC/coursework/Code /Users/tianwang/dropbox/HPC/coursework/Code/build /Users/tianwang/dropbox/HPC/coursework/Code/build /Users/tianwang/dropbox/HPC/coursework/Code/build/CMakeFiles/myprog.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/myprog.dir/depend

