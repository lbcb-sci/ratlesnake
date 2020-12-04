cmake_minimum_required(VERSION 3.2)

project(ratlesnake VERSION 0.4.0
                   LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if (NOT TARGET bioparser)
  add_subdirectory(vendor/ram/vendor/bioparser EXCLUDE_FROM_ALL)
endif ()
if (NOT TARGET ram)
  add_subdirectory(vendor/ram EXCLUDE_FROM_ALL)
endif ()

add_executable(${PROJECT_NAME}
  src/main.cpp)
target_link_libraries(${PROJECT_NAME}
  bioparser
  ram)
target_compile_definitions(${PROJECT_NAME} PRIVATE
  RATLESNAKE_VERSION="v${PROJECT_VERSION}")

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_BINDIR})
