cmake_minimum_required(VERSION 3.2)
project(ratlesnake LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT TARGET bioparser)
    add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
endif()

add_executable(ratlesnake
    src/main.cpp)

target_link_libraries(ratlesnake bioparser)

install(TARGETS ratlesnake DESTINATION ${CMAKE_INSTALL_LIBDIR})
