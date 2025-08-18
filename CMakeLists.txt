cmake_minimum_required(VERSION 3.19)
project(dode LANGUAGES CXX)

# Use C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

## -------------------------- USE MFEM -----------------------------------------
## MFEM is used for the finite element discretization, computation of mass and 
## stiffness matrix, and general sparse matrix-vector operations. The cmake 
## code of this section was provided by MFEM / LLNL directly.

# Use MFEM in this project.
# Import MFEM. The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
#  * mfem_DIR - absolute path to where MFEMConfig.cmake is.
message(STATUS "** MFEM **")
message(STATUS " -> Looking for mfem ...")

# Set local MFEM path (no ~, use full path or $ENV{HOME})
set(MFEM_DIR "mfem-4.8/build")

# Find MFEM in local path only
find_package(MFEM REQUIRED PATHS "${MFEM_DIR}" "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)

message(STATUS " -> Found mfem config in : ${MFEM_DIR} (version ${MFEM_VERSION})")

# Use the same C++ compiler as MFEM
message(STATUS " -> Mfem compiler is     : ${MFEM_CXX_COMPILER}")
if (NOT CMAKE_CXX_COMPILER AND MFEM_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER "${MFEM_CXX_COMPILER}")
endif()

message(STATUS " -> Mfem include dir is  : ${MFEM_INCLUDE_DIRS}")
include_directories(${MFEM_INCLUDE_DIRS})
message(STATUS " -> Mfem library is      : ${MFEM_LIBRARIES}")

# Copy the data directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

## ---------------------------  Compile Static library  ------------------------

# Compile lib${PROJECT_NAME}.a for the solvers defined in src/*. This is a static
# library that can be linked to other executables (e.g. unit-tests and applications).
include_directories(src)
# Get all header and source files in src (and subdirectories)
file(GLOB_RECURSE PROJECT_SOLVERS_CPP src/*.cpp)
file(GLOB_RECURSE PROJECT_SOLVERS_H src/*.hpp)

message(STATUS "** Compile Static Library **")
message(STATUS " -> Adding static library: libdode")
foreach(PROJECT_SOLVER_CPP ${PROJECT_SOLVERS_CPP})
  message(STATUS "    -> ${PROJECT_SOLVER_CPP}")
endforeach()
foreach(PROJECT_SOLVER_H ${PROJECT_SOLVERS_H})
  message(STATUS "    -> ${PROJECT_SOLVER_H}")
endforeach()

add_library(dode STATIC ${PROJECT_SOLVERS_CPP} ${PROJECT_SOLVERS_H})
target_link_libraries(dode PRIVATE ${MFEM_LIBRARIES})

## ---------------------------  Compile Applications  ------------------------
# Compile an executable for each cpp file in applications/*.cpp. Link the
# executable with the static library defined above and with MFEM.
message(STATUS "** Compile Targets **")
file(GLOB PROJECT_SOURCES applications/*.cpp)
foreach(PROJECT_SOURCE ${PROJECT_SOURCES})
  get_filename_component(PROJECT_NAME ${PROJECT_SOURCE} NAME_WE)
  message(STATUS " -> Adding executable  : ${PROJECT_NAME}")
  add_executable(${PROJECT_NAME} ${PROJECT_SOURCE})
  target_link_libraries(${PROJECT_NAME} PRIVATE dode ${MFEM_LIBRARIES})
endforeach()
