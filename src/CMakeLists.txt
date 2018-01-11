# CMake project file for EIRENE

cmake_minimum_required(VERSION 2.8)
project(EIRENE Fortran)

find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})

# Adds compiler flags to all build types.
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "-real-size 64 -fp-model source")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -warn all,nointerfaces -check all")
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-fdefault-real-8 -fdefault-double-8 -fno-range-check")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -fcheck=all")
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
  set(CMAKE_Fortran_FLAGS "-r8")
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

message(STATUS   " --  Setting interface and user routines:       ")
message(STATUS   "       if flag -Dinterface is unspecified, dummies will")
message(STATUS   "       be set. Try -Dinterface=ITER for ITER options.  ")

if(interface)
  get_filename_component(interface ${interface} NAME)
  message(STATUS " --  interface and user routines set to " ${interface})
else()
  message(STATUS " --  Flag -Dinterface unset; will use dummy routines")
endif()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


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

# Decide whether to use plotting or plot_dummy
set(plotpath "")
if( LibGRS AND LibGKS AND LibX11 AND Libz AND LibXt )
  list(APPEND plotpath ${src_plotting})
else()
  list(APPEND plotpath ${src_plot_dummy})
endif()

# Decide whether to use ITER interface or dummy placeholders
# (easily expandable for further interfaces)
set(couplepath "")
if( interface )
  list(APPEND couplepath ${src_interfaces/couple_SOLPS-ITER})
else()
  list(APPEND couplepath ${src_interfaces/couple_Dummy})
endif()
set(userpath "")
if( interface )
  list(APPEND userpath ${src_user-routines/user_iter})
else()
  list(APPEND userpath ${src_user-routines/user_default})
endif()

# Adds EIRENE library to the project.
add_library(EIRENE ${src_assistant} ${src_broadcast} ${src_diagno} 
                   ${src_file-handling} ${src_geometry} 
                   ${src_geometry/time-routines} ${couplepath} 
                   ${src_iterate} ${src_main-routines} ${src_mathematics} 
                   ${src_modules} ${src_output} ${src_particle-tracing} 
                   ${src_photons} ${plotpath} ${src_sampling} 
                   ${src_scoring} ${src_startup-routines} 
                   ${src_surface-processes} ${src_tetrahedra} 
                   ${userpath} ${src_volume-processes} 
                   ${src_xdr-dummy})

# Adds eirene executable to the project and links it to the MPI library.
add_executable(eirene main-routines/eirene_main.f)
target_link_libraries(eirene EIRENE ${MPI_Fortran_LIBRARIES} ${LibGRS} 
					${LibGKS} ${LibX11} ${Libz} ${LibXt})

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

