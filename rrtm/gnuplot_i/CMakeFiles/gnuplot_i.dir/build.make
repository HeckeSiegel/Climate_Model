# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm"

# Include any dependencies generated for this target.
include gnuplot_i/CMakeFiles/gnuplot_i.dir/depend.make

# Include the progress variables for this target.
include gnuplot_i/CMakeFiles/gnuplot_i.dir/progress.make

# Include the compile flags for this target's objects.
include gnuplot_i/CMakeFiles/gnuplot_i.dir/flags.make

gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o: gnuplot_i/CMakeFiles/gnuplot_i.dir/flags.make
gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o: gnuplot_i/gnuplot_i.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o"
	cd "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/gnuplot_i" && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o   -c "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/gnuplot_i/gnuplot_i.c"

gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/gnuplot_i.dir/gnuplot_i.c.i"
	cd "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/gnuplot_i" && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/gnuplot_i/gnuplot_i.c" > CMakeFiles/gnuplot_i.dir/gnuplot_i.c.i

gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/gnuplot_i.dir/gnuplot_i.c.s"
	cd "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/gnuplot_i" && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/gnuplot_i/gnuplot_i.c" -o CMakeFiles/gnuplot_i.dir/gnuplot_i.c.s

gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o.requires:

.PHONY : gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o.requires

gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o.provides: gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o.requires
	$(MAKE) -f gnuplot_i/CMakeFiles/gnuplot_i.dir/build.make gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o.provides.build
.PHONY : gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o.provides

gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o.provides.build: gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o


# Object files for target gnuplot_i
gnuplot_i_OBJECTS = \
"CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o"

# External object files for target gnuplot_i
gnuplot_i_EXTERNAL_OBJECTS =

gnuplot_i/libgnuplot_i.a: gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o
gnuplot_i/libgnuplot_i.a: gnuplot_i/CMakeFiles/gnuplot_i.dir/build.make
gnuplot_i/libgnuplot_i.a: gnuplot_i/CMakeFiles/gnuplot_i.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library libgnuplot_i.a"
	cd "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/gnuplot_i" && $(CMAKE_COMMAND) -P CMakeFiles/gnuplot_i.dir/cmake_clean_target.cmake
	cd "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/gnuplot_i" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gnuplot_i.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
gnuplot_i/CMakeFiles/gnuplot_i.dir/build: gnuplot_i/libgnuplot_i.a

.PHONY : gnuplot_i/CMakeFiles/gnuplot_i.dir/build

gnuplot_i/CMakeFiles/gnuplot_i.dir/requires: gnuplot_i/CMakeFiles/gnuplot_i.dir/gnuplot_i.c.o.requires

.PHONY : gnuplot_i/CMakeFiles/gnuplot_i.dir/requires

gnuplot_i/CMakeFiles/gnuplot_i.dir/clean:
	cd "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/gnuplot_i" && $(CMAKE_COMMAND) -P CMakeFiles/gnuplot_i.dir/cmake_clean.cmake
.PHONY : gnuplot_i/CMakeFiles/gnuplot_i.dir/clean

gnuplot_i/CMakeFiles/gnuplot_i.dir/depend:
	cd "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm" "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/gnuplot_i" "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm" "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/gnuplot_i" "/home/s/S.Legler/Documents/Advanced Atmospheric Physics/Climate_Model/rrtm/gnuplot_i/CMakeFiles/gnuplot_i.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : gnuplot_i/CMakeFiles/gnuplot_i.dir/depend

