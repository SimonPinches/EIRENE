# CMake project file for EIRENE

cmake_minimum_required(VERSION 2.8)

project(EIRENE Fortran)

# Set a default build type if none was specified.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Release" "Debug" "Coverage")
  # "Coverage" only available for GNU compiler
endif()

# Adds compiler flags to all build types.
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS_RELEASE "-fp-model source -fpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -fpp -gen-interfaces -warn all,interfaces -check all -traceback")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -cpp -fopenmp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -cpp -fopenmp -Wall -fcheck=all -Wdo-subscript -Wtarget-lifetime")
  set(CMAKE_Fortran_FLAGS_COVERAGE "-g -fprofile-arcs -ftest-coverage")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -cpp -mp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -cpp -mp -C -Mchkptr -Mchkstk")
endif()

message("-- Adding Fortran 2003 flag...")
add_definitions(-DF2003)



# Output names and postfix.
set(EIRENE_LIBRARY_NAME "EIRENE" CACHE STRING "EIRENE library name")
set(EIRENE_EXECUTABLE_NAME "eirene" CACHE STRING "EIRENE executable name")
set(EIRENE_POSTFIX "" CACHE STRING "EIRENE output names postfix")

# Check for modules and library for JSON
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (NOT JSON_MODULES)
  message(STATUS "Setting default directory for JSON_MODULES.")
  set(JSON_MODULES "/var/lib")
endif()
if(JSON_MODULES)
  get_filename_component(JSON_MODULES ${JSON_MODULES} ABSOLUTE)
  include_directories(${JSON_MODULES})
  message(STATUS "  JSON_MODULES:   " ${JSON_MODULES})
else()
#  message(FATAL_ERROR "Directory of JSON module files not set. "
#                      "Use cmake -DJSON_MODULES=\"<directory>\" to set it.")
  message(STATUS "Directory of JSON module files not set. "
                      "Use cmake -DJSON_MODULES=\"<directory>\" to set it.")
endif()
if (NOT LibJSON)
  message(STATUS "Setting default directory for LibJSON.")
  set(LibJSON "/var/lib")
endif()

if(LibJSON)
  get_filename_component(LibJSON ${LibJSON} ABSOLUTE)
  message(STATUS "  Lib json_fortran:   " ${LibJSON})
else()
  message(STATUS  "Path to LibJSON unset")
  message(STATUS  "You may want to try -DLibJSON=$PATH_TO_JSON_LIB/libjsonfortran.so")
endif()

# Libraries used for plotting
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
message(STATUS    " --  Setting libraries for internal plotting:  --  --")

if(LibGKS)
  get_filename_component(LibGKS ${LibGKS} ABSOLUTE)
else()
  message(STATUS  "Path to LibGKS unset")
  message(STATUS  "You may want to try -DLibGKS=$PATH_TO_GKS_LIB/libgks.a")
endif()

if(LibGRS)
  get_filename_component(LibGRS ${LibGRS} ABSOLUTE)
else()
  message(STATUS  "Path to LibGRS unset")
  message(STATUS  "You may want to try -DLibGRS=$PATH_TO_GRSOFT_LIB/libgr.a")
endif()

if(LibX11)
  get_filename_component(LibX11 ${LibX11} ABSOLUTE)
else()
  message(STATUS  "Path to LibX11 unset")
  message(STATUS  "You may want to try -DLibX11=$PATH_TO_X11_LIB/libX11.so")
endif()

if(Libz)
  get_filename_component(Libz ${Libz} ABSOLUTE)
else()
  message(STATUS  "Path to Libz unset")
  message(STATUS  "You may want to try -DLibz=$PATH_TO_z_LIB/libz.so")
endif()

if(LibXt)
  get_filename_component(LibXt ${LibXt} ABSOLUTE)
else()
  message(STATUS  "Path to LibXt unset")
  message(STATUS  "You may want to try -DLibXt=$PATH_TO_Xt_LIB/libXt.so")
endif()

if(LibGKS AND LibGRS AND LibX11 AND Libz AND LibXt)
  message(STATUS "Found the following libraries used for plotting:")
  message(STATUS "  Lib GKS: " ${LibGKS})
  message(STATUS "  Lib GRS: " ${LibGRS})
  message(STATUS "  Lib X11: " ${LibX11})
  message(STATUS "  Lib Xt:  " ${LibXt})
  message(STATUS "  Lib z:   " ${Libz})
  
  message(STATUS "For running, please export the following 2 statements "
		 "into your console:")
  message(STATUS "export GLI_HOME=$PATH_TO_GKS_LIB")
  message(STATUS "export GRSOFT_DEVICE=62")
  message(STATUS " --  Plotting will be enabled.")
else()
  message(STATUS " --  Plotting will be disabled.")
endif()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Control of interface and user-routines to compile.
set(EIRENE_INTERFACE "Dummy" CACHE STRING "EIRENE interface to compile")
set(EIRENE_USER-ROUTINES "default" CACHE STRING "EIRENE user-routines to compile")

message(STATUS   " --  Setting interface and user-routines to compile:")
message(STATUS   "       If flag -DEIRENE_INTERFACE is unspecified, \"Dummy\"")
message(STATUS   "       will be set. Try e.g. -DEIRENE_INTERFACE=SOLPS-ITER.  ")
message(STATUS   "       If flag -DEIRENE_USER-ROUTINES is unspecified, \"default\"")
message(STATUS   "       will be set. Try e.g. -DEIRENE_USER-ROUTINES=iter.  ")

# Add option to control whether trace output will be compiled.
option(TRACE "Build with trace output" OFF)
if (TRACE)
  message("-- Adding TRACE flag...")
  add_definitions(-DTRACE)
endif()

# Add option to control whether binary files are also written in ASCII.
option(CHECKBIN "Check binary files" OFF)
if (CHECKBIN)
  message("-- Adding CHECKBIN flag...")
  add_definitions(-DCHECKBIN)
endif()

# Add option to control whether MPI will be used.
option(MPI "Build with MPI" ON)
if (MPI)
  message("-- Adding USE_MPI flag...")
  add_definitions(-DUSE_MPI)
  find_package(MPI REQUIRED)
#  add definitions(-DMPI_VERSION=${MPI_VERSION})
  include_directories(${JSON_MODULES} ${MPI_Fortran_INCLUDE_PATH})
endif()

# Add option to enable OpenMP with the parallel region defined in eirmod_mcarlo.
option(OPENMP "Build for use with internal parallel region" OFF)
if (OPENMP)
  cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
  message("-- Adding USE_OPENMP flag...")
  add_definitions(-DUSE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# Add option to enable OpenMP with the parallel region defined in the calling program.
option(EXT_OPENMP "Build for use with external parallel region" OFF)
if (EXT_OPENMP)
  cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
  message("-- Adding USE_EXT_OPENMP flag...")
  add_definitions(-DUSE_OPENMP)
  add_definitions(-DUSE_EXT_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# Adds openMP using the old style CMake approach to ensure backwards compatibility
# This should be changed at some point to use the newer OpenMP Fortran capabilities
if (OPENMP_FOUND)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

# Finds all fortran files (*.[fF]/*.[fF]90) in the source directory and extract the 
# relative path. Relative path is stored in variable 'fortran_dirs'.
message(STATUS "Detecting fortran files (*.[fF]/*.[fF]90) in subdirectories of " 
               ${PROJECT_SOURCE_DIR})
file(GLOB_RECURSE fortran_files RELATIVE ${PROJECT_SOURCE_DIR} 
     ${PROJECT_SOURCE_DIR}/*.[fF]
     ${PROJECT_SOURCE_DIR}/*.[fF]90)
set(fortran_dirs "")
foreach(fortran_file ${fortran_files})
  get_filename_component(fortran_dir ${fortran_file} PATH)
  list(APPEND fortran_dirs ${fortran_dir})
endforeach()
list(REMOVE_DUPLICATES fortran_dirs)

# Sets variables 'src_[relative path to subdirectory]' which contain a list of 
# all fortran files (*.[fF]/*.[fF]90) in this subdirectory.
message(STATUS "Configuring variables src_[relative path to subdirectory]")
foreach(fortran_dir ${fortran_dirs})
  file(GLOB src_${fortran_dir} ${PROJECT_SOURCE_DIR}/${fortran_dir}/*.[fF]
                               ${PROJECT_SOURCE_DIR}/${fortran_dir}/*.[fF]90)
endforeach()

# Remove program file to be able to add it seperately if required.
list(GET src_main-routines 0 fortran_file)
get_filename_component(fortran_dir ${fortran_file} PATH)
list(REMOVE_ITEM src_main-routines ${fortran_dir}/eirene_main.f)

# Decide whether to use plotting or plot_dummy.
set(plotpath "")
if( LibGRS AND LibGKS AND LibX11 AND Libz AND LibXt )
  list(APPEND plotpath ${src_plotting})
else()
#  list(APPEND plotpath ${src_plot_dummy})
  list(APPEND plotpath ${src_plotting})
  list(APPEND plotpath ${src_grsoft-dummy})
endif()

# Specify interface and user-routines to use.
set(src_interface/x ${src_interfaces/couple_${EIRENE_INTERFACE}})
set(src_user-routines/x ${src_user-routines/user_${EIRENE_USER-ROUTINES}})
if( NOT src_interface/x )
  message(FATAL_ERROR "No interface named \"${EIRENE_INTERFACE}\"")
endif()
if( NOT src_user-routines/x )
  message(FATAL_ERROR "No user-routines named \"${EIRENE_USER-ROUTINES}\"")
endif()

# Adds EIRENE library to the project.
add_library(EIRENE ${src_assistant} ${src_broadcast} ${src_diagno} 
                   ${src_file-handling} ${src_function_modules} ${src_geometry}
                   ${src_geometry/time-routines} ${src_interface/x} 
                   ${src_iterate} ${src_main-routines} ${src_mathematics} 
                   ${src_modules} ${src_output} ${src_particle-tracing} 
                   ${src_photons} ${plotpath} ${src_sampling} 
                   ${src_scoring} ${src_startup-routines}  
                   ${src_surface-processes} ${src_tetrahedra} 
                   ${src_user-routines/x} ${src_volume-processes})

# Adds eirene executable to the project and links it to the MPI library.
add_executable(eirene main-routines/eirene_main.f)
if (MPI)
  target_link_libraries(eirene EIRENE ${LibJSON} ${MPI_Fortran_LIBRARIES})
else()
  target_link_libraries(eirene EIRENE ${LibJSON})
endif()

if(LibGKS AND LibGRS AND LibX11 AND Libz AND LibXt AND LibJSON)
  target_link_libraries(eirene EIRENE ${LibGRS} ${LibGKS} ${LibX11} ${LibXt} ${LibJSON})
endif()

# This is the newer approach for adding OpenMP - ver >= 3.9
#if(OPENMP OR EXT_OPENMP)
#  target_compile_options(EIRENE PRIVATE ${OpenMP_Fortran_FLAGS})
#  target_compile_options(eirene PRIVATE ${OpenMP_Fortran_FLAGS})
#  target_link_libraries(eirene EIRENE ${OpenMP_Fortran_LIBRARIES})
#endif()

# Sets output directory of the library and executbale and output name 
# of the executable.
set_target_properties(EIRENE PROPERTIES
                      ARCHIVE_OUTPUT_DIRECTORY 
                      ${PROJECT_BINARY_DIR}/../lib${CMAKE_BUILD_TYPE}
                      OUTPUT_NAME ${EIRENE_LIBRARY_NAME}${EIRENE_POSTFIX})
set_target_properties(eirene PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY 
                      ${PROJECT_BINARY_DIR}/../bin${CMAKE_BUILD_TYPE}
                      OUTPUT_NAME ${EIRENE_EXECUTABLE_NAME}${EIRENE_POSTFIX})

# Gets the git revision hash and stores it in the variable 'git_hash'.
find_package(Git)

if(GIT_FOUND)
  if(EXISTS ${PROJECT_SOURCE_DIR}/../.git)
    set(git_hash_command log --pretty=format:"%ci" --abbrev-commit -1)
    execute_process(COMMAND ${GIT_EXECUTABLE} ${git_hash_command} 
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/.. 
                    OUTPUT_VARIABLE git_date 
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    set(git_hash_command log --pretty=format:"%H" --abbrev-commit -1)
    execute_process(COMMAND ${GIT_EXECUTABLE} ${git_hash_command} 
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/.. 
                    OUTPUT_VARIABLE git_hash 
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    set(git_hash_command log --pretty=format:"%d" --abbrev-commit -1)
    execute_process(COMMAND ${GIT_EXECUTABLE} ${git_hash_command} 
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/.. 
                    OUTPUT_VARIABLE git_branch 
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    set(git_hash_command log --pretty=format:"%ci %H %d" --abbrev-commit -1)
    execute_process(COMMAND ${GIT_EXECUTABLE} ${git_hash_command} 
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/.. 
                    OUTPUT_VARIABLE git_full_hash 
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    # Configures the version CMakeLists file.
    configure_file(${PROJECT_SOURCE_DIR}/cmake/version.cmake.in
                   version.cmake
                   @ONLY)

    # Trigger the EIRENE version output.
    add_custom_target(version 
                      ${CMAKE_COMMAND} -D PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                                       -P version.cmake)
    add_dependencies(EIRENE version)
    add_dependencies(eirene version)
  endif()
endif()

# Adds a target to generate a documentation with Doxygen.
find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(${PROJECT_SOURCE_DIR}/cmake/Doxyfile.in 
                 ${PROJECT_BINARY_DIR}/../doc/Doxyfile @ONLY)

  add_custom_target(doc
                   ${DOXYGEN_EXECUTABLE} Doxyfile
                   WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/../doc
                   COMMENT "Generating documentation with Doxygen" VERBATIM)

  message("-----------------------------------------------------")
  message("To generate the documentation, type \"make doc\"")
  message("-----------------------------------------------------")
endif()

