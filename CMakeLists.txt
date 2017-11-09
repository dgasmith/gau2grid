cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(gau2grid
        VERSION 0.1
        LANGUAGES C)
set(gau2grid_AUTHORS      "Daniel G. A. Smith")
set(gau2grid_DESCRIPTION  "Fast computation of a gaussian and its derivative on a grid")
set(gau2grid_URL          "https://github.com/dgasmith/gau2grid")
set(gau2grid_LICENSE      "BSD 3-clause")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#############################  Options: Build How?  #############################
include(psi4OptionsTools)
#option_with_flags(ENABLE_XHOST "Enables processor-specific optimization" ON
#                  "-xHost" "-march=native")
option_with_default(CMAKE_BUILD_TYPE "Build type (Release or Debug)" Release)
#option_with_default(CMAKE_INSTALL_LIBDIR "Directory to which libraries installed" lib)
#option_with_default(PYMOD_INSTALL_LIBDIR "Location within CMAKE_INSTALL_LIBDIR to which python modules are installed" /)
#option_with_default(ENABLE_GENERIC "Enables mostly static linking of system libraries for shared library" OFF)
option_with_default(gau2grid_CXX_STANDARD "Specify C++ standard for core gau2grid" 11)

########################  Process & Validate Options  ##########################
#include(GNUInstallDirs)
include(autocmake_safeguards)
include(custom_color_messages)

#if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#    set(CMAKE_INSTALL_PREFIX "/usr/local/gau2grid" CACHE PATH "Install path" FORCE)
#endif()
#message(STATUS "gau2grid install: ${CMAKE_INSTALL_PREFIX}")

set(CMAKE_CXX_STANDARD "${gau2grid_CXX_STANDARD}")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#  <<  Pybind11 & Python  >>
set(PYBIND11_CPP_STANDARD "-std=c++${CMAKE_CXX_STANDARD}")
find_package(pybind11 2.0.0 CONFIG REQUIRED)
message(STATUS "${Cyan}Using pybind11${ColourReset}: ${pybind11_INCLUDE_DIR} (version ${pybind11_VERSION} for Py${PYTHON_VERSION_STRING} and ${PYBIND11_CPP_STANDARD})")
message(STATUS "${Cyan}Using Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}${ColourReset}: ${PYTHON_EXECUTABLE}")

## external projects manage their own OpenMP and c++YY flags, so only add to CXX_FLAGS for psi4-core
#include(autocmake_omp)
#include(custom_cxxstandard)
#include(custom_static_library)

################################  Main Project  ################################

add_custom_command(
    OUTPUT  gau2grid_phi.c gau2grid_phi_grad.c
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys; \
                                     sys.path.append('${CMAKE_SOURCE_DIR}'); \
                                     import gau2grid as gg; \
                                     gg.c_gen.generate_c_gau2grid(3, path='${CMAKE_CURRENT_BINARY_DIR}')"
    DEPENDS gau2grid/c_generator.py
    VERBATIM)

set(sources_list ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_phi.c)
                 #main.cpp)

add_executable(hello ${sources_list})
add_library(gau2grid ${sources_list})

pybind11_add_module(pygau2grid NO_EXTRAS ${sources_list})
target_link_libraries(pygau2grid PRIVATE gau2grid)
target_link_libraries(pygau2grid PUBLIC pybind11::module)

