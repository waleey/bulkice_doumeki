#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21 FATAL_ERROR)
project(bulkice_doumeki)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#set(BOOST_ROOT "~/local/boost_1_xx_0")
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
	    #add_definitions(-DG4VIS_USE_RAYTRACERX)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
#Copying the data program to the build folder
#
# Set the source directory you want to copy
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/../data")

# Set the destination directory where you want to copy the source directory
set(DESTINATION_DIR "${PROJECT_SOURCE_DIR}/build/")

# Copy the directory
file(COPY "${SOURCE_DIR}" DESTINATION "${DESTINATION_DIR}")


#----------------------------------------------------------------------------
#Copying the merger program to the build folder
#
# Set the source directory you want to copy
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/../merger")

# Set the destination directory where you want to copy the source directory
set(DESTINATION_DIR "${PROJECT_SOURCE_DIR}/build/")

# Copy the directory
file(COPY "${SOURCE_DIR}" DESTINATION "${DESTINATION_DIR}")

#----------------------------------------------------------------------------
#Creating the output folder
#
#set(OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/build/output")

# Add a custom command to create the folder
#add_custom_command(
#   OUTPUT "${OUTPUT_DIRECTORY}"
#   COMMAND ${CMAKE_COMMAND} -E make_directory "${OUTPUT_DIRECTORY}"
#   COMMENT "Creating directory ${OUTPUT_DIRECTORY}"
#)

# Set a target that depends on the custom command
#add_custom_target(create_my_folder DEPENDS "${OUTPUT_DIRECTORY}")

# Add dependencies to make sure the folder is created before building
#add_dependencies(${PROJECT_NAME} create_my_folder)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
#set(OProfile_DIR "/bin//operf")
# Link against OProfile libraries
#find_package(OProfile REQUIRED)
#include_directories(${OPROFILE_INCLUDE_DIRS})
#target_link_libraries(bulkice_doumeki ${OPROFILE_LIBRARIES})
add_executable(bulkice_doumeki bulkice_doumeki.cc ${sources} ${headers} ${TOOLS_FORTRAN_OBJECTS})
target_link_libraries(bulkice_doumeki ${Geant4_LIBRARIES} ${HBOOK_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build TestEm1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(bulkice_doumeki_SCRIPTS
    100run_pos.mac vis.mac run_pos.mac ray_trace.mac 10run_pos.mac env.sh 
  )
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
# Link against OProfile libraries
#find_package(OProfile REQUIRED)
#include_directories(${OPROFILE_INCLUDE_DIRS})
#target_link_libraries(your_target_name ${OPROFILE_LIBRARIES})

foreach(_script ${bulkice_doumeki_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS bulkice_doumeki DESTINATION bin)

