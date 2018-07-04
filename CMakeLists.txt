cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(gau2grid
        VERSION 1.1.0
        LANGUAGES C)
set(gau2grid_AUTHORS      "Daniel G. A. Smith")
set(gau2grid_DESCRIPTION  "Fast computation of a gaussian and its derivative on a grid")
set(gau2grid_URL          "https://github.com/dgasmith/gau2grid")
set(gau2grid_LICENSE      "BSD 3-clause")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


#############################  Options: Build How?  #############################
include(psi4OptionsTools)
option_with_default(MAX_AM "The maximum gaussian angular momentum to compile" 8)
option_with_default(SPHERICAL_ORDER "The order for the spherical basis functions. Options: cca, gaussian" "gaussian")
option_with_default(CARTESIAN_ORDER "The order for the cartesian basis functions. Options: row" "row")

option_with_default(MAX_AM "The maximum gaussian angular momentum to compile" 8)
option_with_default(CMAKE_BUILD_TYPE "Build type (Release or Debug)" Release)
option_with_flags(ENABLE_XHOST "Enables processor-specific optimization" ON
                  "-xHost" "-march=native")

option_with_default(BUILD_FPIC "Libraries will be compiled with position independent code" ON)
option_with_print(BUILD_SHARED_LIBS "Build final library as shared, not static" ON)
option_with_default(ENABLE_GENERIC "Enables mostly static linking of system libraries for shared library" OFF)

# Warnings
if((${BUILD_SHARED_LIBS}) AND NOT ${BUILD_FPIC})
    message(FATAL_ERROR "BUILD_SHARED_LIBS ON and BUILD_FPIC OFF are incompatible, as shared library requires position independent code")
endif()
if(NOT ((${SPHERICAL_ORDER} STREQUAL "cca") OR (${SPHERICAL_ORDER} STREQUAL "gaussian")))
    message(FATAL_ERROR "SPHERICAL_ORDER must be 'cca' or 'gaussian': ${SPHERICAL_ORDER}")
endif()

# Install
option_with_default(CMAKE_INSTALL_LIBDIR "Directory to which libraries installed" lib)
option_with_default(PYMOD_INSTALL_LIBDIR "Location within CMAKE_INSTALL_LIBDIR to which python modules are installed
                                          Must start with: / . Used to imitate python install: /python3.6/site-packages ." /)
option_with_print(INSTALL_PYMOD "Additionally installs as independent python module in PYMOD_INSTALL_LIBDIR" OFF)

########################  Process & Validate Options  ##########################
include(autocmake_safeguards)
include(custom_color_messages)
include(custom_static_library)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/usr/local/gau2grid" CACHE PATH "Install path" FORCE)
endif()
message(STATUS "gau2grid install: ${CMAKE_INSTALL_PREFIX}")

#  <<  Python  >>
set(Python_ADDITIONAL_VERSIONS 3.7 3.6 3.5)  # adjust with CMake minimum FindPythonInterp
find_package(PythonLibsNew 2.7 REQUIRED)
message(STATUS "${Cyan}Found Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}${ColourReset}: ${PYTHON_EXECUTABLE} (found version ${PYTHON_VERSION_STRING})")


################################  Main Project  ################################
add_custom_command(
    OUTPUT  gau2grid.h gau2grid_phi.c gau2grid_deriv1.c gau2grid_deriv2.c gau2grid_spherical.c gau2grid_helper.c
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys; \
                                     sys.path.append('${CMAKE_SOURCE_DIR}'); \
                                     import gau2grid as gg; \
                                     gg.c_gen.generate_c_gau2grid(${MAX_AM}, path='${CMAKE_CURRENT_BINARY_DIR}', cartesian_order='${CARTESIAN_ORDER}', spherical_order='${SPHERICAL_ORDER}')"
    DEPENDS gau2grid/c_generator.py
    VERBATIM)

set(sources_list ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_phi.c
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_deriv1.c
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_deriv2.c
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_spherical.c
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_helper.c)

add_library(gg ${sources_list})
set_target_properties(gg PROPERTIES COMPILE_FLAGS "-std=c99")
set_target_properties(gg PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC})

if(${BUILD_SHARED_LIBS})
    target_link_libraries(gg PRIVATE ${LIBC_INTERJECT})
endif()


###################################  Install  ##################################
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(PN ${PROJECT_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gau2grid.h
              ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_pragma.h
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

# Install our files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(EXPORT "${PN}Targets"
        NAMESPACE "${PN}::"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

if(${INSTALL_PYMOD})
    install(DIRECTORY gau2grid
            DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}
            USE_SOURCE_PERMISSIONS
            FILES_MATCHING PATTERN "*.py" PATTERN "*.npz")

    install(FILES $<TARGET_FILE:gg>
            DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/gau2grid)
endif()
