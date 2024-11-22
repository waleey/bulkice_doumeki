#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21 FATAL_ERROR)
project(bulkice_doumeki)

#----------------------------------------------------------------------------
# Find Geant4 package with UI and Visualization
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Boost integration
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})
add_definitions(-DG4VIS_USE_OPENGL_X11 -DG4VIS_USE_VERBOSE)

#----------------------------------------------------------------------------
# Locate sources and headers
include_directories(${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Explicit OpenGL and X11 Linking (Optional)
find_package(OpenGL REQUIRED)
find_package(X11 REQUIRED)
include_directories(${X11_INCLUDE_DIR})
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Copy resources and create output directories
file(COPY "${PROJECT_SOURCE_DIR}/../data" DESTINATION "${PROJECT_SOURCE_DIR}/build/")
file(COPY "${PROJECT_SOURCE_DIR}/../merger" DESTINATION "${PROJECT_SOURCE_DIR}/build/")
file(MAKE_DIRECTORY "${PROJECT_SOURCE_DIR}/build/output")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
add_executable(bulkice_doumeki bulkice_doumeki.cc ${sources})
target_link_libraries(bulkice_doumeki ${Geant4_LIBRARIES} ${OPENGL_LIBRARIES} ${X11_LIBRARIES} ${Boost_LIBRARIES})

#----------------------------------------------------------------------------
# Copy scripts to build directory
set(bulkice_doumeki_SCRIPTS
    100run_pos.mac vis.mac run_pos.mac ray_trace.mac 10run_pos.mac env.sh 
)
foreach(_script ${bulkice_doumeki_SCRIPTS})
  configure_file(${PROJECT_SOURCE_DIR}/${_script} ${PROJECT_BINARY_DIR}/${_script} COPYONLY)
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
install(TARGETS bulkice_doumeki DESTINATION bin)