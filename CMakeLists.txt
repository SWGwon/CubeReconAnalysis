cmake_minimum_required(VERSION 3.0)

# Configure the dependencies
find_package(ROOT REQUIRED)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

if(ROOT_Eve_LIBRARY)
  message("Building the cube-analysis event analysislay")
else(ROOT_Eve_LIBRARY)
  message("NOT BUILDING THE EVENT DISPLAY: eve not found")
  return()
endif(ROOT_Eve_LIBRARY)

# Make sure the current directories are available for the root
include_directories(${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../io
  ${CMAKE_CURRENT_SOURCE_DIR}/../src
  ${CMAKE_CURRENT_SOURCE_DIR}/../tools)

set(source 
    src/analysisExample.cxx)

add_executable(analysisExample ${source})
target_link_libraries(analysisExample LINK_PUBLIC cuberecon_io cuberecon_tools ${ROOT_LIBRARIES})
#target_link_libraries(analysisExample)
install(TARGETS analysisExample RUNTIME DESTINATION bin)

# If this is ROOT6 or later, then install the rootmap and pcm files.
if(${ROOT_VERSION} VERSION_GREATER 6)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/liblib_cubeanalysis.g.rootmap
    DESTINATION lib)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/liblib_cubeanalysis.g_rdict.pcm
    DESTINATION lib)
endif(${ROOT_VERSION} VERSION_GREATER 6)
