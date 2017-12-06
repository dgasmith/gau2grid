cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(gau2grid
        VERSION 0.1
        LANGUAGES C CXX)
set(gau2grid_AUTHORS      "Daniel G. A. Smith")
set(gau2grid_DESCRIPTION  "Fast computation of a gaussian and its derivative on a grid")
set(gau2grid_URL          "https://github.com/dgasmith/gau2grid")
set(gau2grid_LICENSE      "BSD 3-clause")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#############################  Options: Build How?  #############################
include(psi4OptionsTools)
option_with_flags(ENABLE_XHOST "Enables processor-specific optimization" ON
                  "-xHost" "-march=native")
option_with_default(CMAKE_BUILD_TYPE "Build type (Release or Debug)" Release)
option_with_print(BUILD_SHARED_LIBS "Build final library as shared, not static" OFF)
option_with_default(BUILD_FPIC "Libraries will be compiled with position independent code" ON)
if((${BUILD_SHARED_LIBS}) AND NOT ${BUILD_FPIC})
    message(FATAL_ERROR "BUILD_SHARED_LIBS ON and BUILD_FPIC OFF are incompatible, as shared library requires position independent code")
endif()
#option_with_default(CMAKE_INSTALL_LIBDIR "Directory to which libraries installed" lib)
#option_with_default(PYMOD_INSTALL_LIBDIR "Location within CMAKE_INSTALL_LIBDIR to which python modules are installed" /)
option_with_default(ENABLE_GENERIC "Enables mostly static linking of system libraries for shared library" OFF)
option_with_default(gau2grid_CXX_STANDARD "Specify C++ standard for core gau2grid" 11)
option_with_default(ENABLE_PYTHON_BINDINGS "Enables pybind11 interface to C++ library" OFF)
option_with_default(MAX_AM_ERI "The maximum gaussian angular momentum to compile" 8)

########################  Process & Validate Options  ##########################
include(autocmake_safeguards)
include(custom_color_messages)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/usr/local/gau2grid" CACHE PATH "Install path" FORCE)
endif()
message(STATUS "gau2grid install: ${CMAKE_INSTALL_PREFIX}")

set(CMAKE_CXX_STANDARD "${gau2grid_CXX_STANDARD}")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#  <<  Pybind11 & Python  >>
set(PYBIND11_CPP_STANDARD "-std=c++${CMAKE_CXX_STANDARD}")
find_package(pybind11 2.0.0 CONFIG REQUIRED)
message(STATUS "${Cyan}Using pybind11${ColourReset}: ${pybind11_INCLUDE_DIR} (version ${pybind11_VERSION} for Py${PYTHON_VERSION_STRING} and ${PYBIND11_CPP_STANDARD})")
message(STATUS "${Cyan}Using Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}${ColourReset}: ${PYTHON_EXECUTABLE}")

#include(custom_cxxstandard)
include(custom_static_library)

################################  Main Project  ################################

add_custom_command(
    OUTPUT  gau2grid_phi.cc gau2grid_deriv1.cc gau2grid_deriv2.cc gau2grid_spherical.cc pygg_core.cc gau2grid_helper.cc
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys; \
                                     sys.path.append('${CMAKE_SOURCE_DIR}'); \
                                     import gau2grid as gg; \
                                     gg.c_gen.generate_c_gau2grid(${MAX_AM_ERI}, path='${CMAKE_CURRENT_BINARY_DIR}')"
    DEPENDS gau2grid/c_generator.py
    VERBATIM)

set(sources_list ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_phi.cc
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_deriv1.cc
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_deriv2.cc
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_spherical.cc
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_helper.cc)

add_library(gg ${sources_list})
set_target_properties(gg PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC})

if(${BUILD_SHARED_LIBS})
    target_link_libraries(gg PRIVATE ${LIBC_INTERJECT})
endif()

if(${ENABLE_PYTHON_BINDINGS})
    #list(APPEND sources_list ${CMAKE_CURRENT_BINARY_DIR}/pygg_core.cc)
    pybind11_add_module(pygg_core NO_EXTRAS ${CMAKE_CURRENT_BINARY_DIR}/pygg_core.cc)
    target_link_libraries(pygg_core PRIVATE gg)
    target_link_libraries(pygg_core PUBLIC pybind11::module)
endif()

###################################  Install  ##################################
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(PN ${PROJECT_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gau2grid.h
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PN})

install(TARGETS gg
        EXPORT "${PN}Targets"
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

target_include_directories(gg INTERFACE
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# GNUInstallDirs "DATADIR" wrong here; CMake search path wants "share".
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}")
configure_package_config_file(cmake/${PN}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
                                 VERSION ${${PN}_VERSION}
                                 COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(EXPORT "${PN}Targets"
        NAMESPACE "${PN}::"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

