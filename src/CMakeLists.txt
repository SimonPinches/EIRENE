# CMake project file for EIRENE

cmake_minimum_required(VERSION 2.8)
project(EIRENE Fortran)

find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})

# Adds compiler flags to all build types.
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "-real-size 64 -fp-model source")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -warn all -check all")
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-fdefault-real-8 -fdefault-double-8 -fno-range-check")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -fcheck=all")
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
  set(CMAKE_Fortran_FLAGS "-r8")
endif()

# Add option to control whether trace output will be compiled
option(TRACE "Build with trace output" OFF)
if (TRACE)
  message("-- Adding TRACE flag...")
  add_definitions(-DTRACE)
endif()

# Finds all fortran files (*.f) in the source directory and extract the 
# relative path. Relative path is stored in variable 'fortran_dirs'.
message(STATUS "Detecting fortran files (*.[fF]) in subdirectories of " 
               ${PROJECT_SOURCE_DIR})
file(GLOB_RECURSE fortran_files RELATIVE ${PROJECT_SOURCE_DIR} 
     ${PROJECT_SOURCE_DIR}/*.[fF])
set(fortran_dirs "")
foreach(fortran_file ${fortran_files})
  get_filename_component(fortran_dir ${fortran_file} PATH)
  list(APPEND fortran_dirs ${fortran_dir})
endforeach()
list(REMOVE_DUPLICATES fortran_dirs)

# Sets variables 'src_[relative path to subdirectory]' which contain a list of 
# all fortran files (*.[fF]) in this subdirectory.
message(STATUS "Configuring variables src_[relative path to subdirectory]")
foreach(fortran_dir ${fortran_dirs})
  file(GLOB src_${fortran_dir} ${PROJECT_SOURCE_DIR}/${fortran_dir}/*.[fF])
endforeach()

# Remove program file to be able to add it seperately if required.
list(GET src_main-routines 0 fortran_file)
get_filename_component(fortran_dir ${fortran_file} PATH)
list(REMOVE_ITEM src_main-routines ${fortran_dir}/eirene_main.f)

# Adds EIRENE library to the project.
add_library(EIRENE ${src_assistant} ${src_broadcast} ${src_diagno} 
                   ${src_file-handling} ${src_geometry} 
                   ${src_geometry/time-routines} ${src_interfaces/couple_Dummy} 
                   ${src_iterate} ${src_main-routines} ${src_mathematics} 
                   ${src_modules} ${src_output} ${src_particle-tracing} 
                   ${src_photons} ${src_plot_dummy} ${src_sampling} 
                   ${src_scoring} ${src_startup-routines} 
                   ${src_surface-processes} ${src_tetrahedrons} 
                   ${src_user-routines/user_default} ${src_volume-processes} 
                   ${src_xdr-dummy})

# Adds eirene executable to the project and links it to the MPI library.
add_executable(eirene main-routines/eirene_main.f)
target_link_libraries(eirene EIRENE ${MPI_Fortran_LIBRARIES})

# Sets output directory of the library and executbales.
set_target_properties(EIRENE PROPERTIES
                      ARCHIVE_OUTPUT_DIRECTORY 
                      ${PROJECT_BINARY_DIR}/../lib${CMAKE_BUILD_TYPE})
set_target_properties(eirene PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY 
                      ${PROJECT_BINARY_DIR}/../bin${CMAKE_BUILD_TYPE})

# Gets the git revision hash and stores it in the variable 'git_hash'.
find_package(Git)
if(GIT_FOUND)
  if(EXISTS ${PROJECT_SOURCE_DIR}/../.git)
    set(git_hash_command log --pretty=format:"%ci %H %d" --abbrev-commit -1)
    execute_process(COMMAND ${GIT_EXECUTABLE} ${git_hash_command} 
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/.. 
                    OUTPUT_VARIABLE git_hash 
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
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
  message("To generate the documentation, type \"make docs\"")
  message("-----------------------------------------------------")
endif()

