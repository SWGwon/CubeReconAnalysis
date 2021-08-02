cmake_minimum_required(VERSION 3.0)

# Configure the dependencies
find_package(ROOT REQUIRED)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

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
install(TARGETS analysisExample RUNTIME DESTINATION bin)

find_package(EDepSim)
add_executable(EventSelection src/EventSelection.cxx src/EventSelection.hxx)
target_link_libraries(EventSelection LINK_PUBLIC cuberecon_io cuberecon_tools ${ROOT_LIBRARIES} EDepSim::edepsim_io)
install(TARGETS EventSelection RUNTIME DESTINATION bin)
