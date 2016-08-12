cmake_minimum_required(VERSION 3.0)
project(libefp LANGUAGES C)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

################################### Options ####################################
include(Psi4Macros)
option_with_default(CMAKE_BUILD_TYPE "Build type" Release)
option_with_print(BUILD_SHARED_LIBS "Build final library as shared, not static" OFF)
option_with_default(BUILD_FPIC "Libraries will be compiled with position independent code" ON)
if(${BUILD_SHARED_LIBS} AND NOT ${BUILD_FPIC})
    message(FATAL_ERROR "BUILD_SHARED_LIBS ON and BUILD_FPIC OFF are incompatible, as shared library requires position independent code")
endif()
option_with_print(ENABLE_OPENMP "Enable OpenMP parallelization" ON)
option_with_print(ENABLE_GENERIC "Enable mostly static linking in shared library" OFF)
option_with_flags(ENABLE_XHOST "Enable processor-specific optimization" ON
                  "-xHost" "-march=native")
option_with_print(FRAGLIB_UNDERSCORE_L "Installed fragment library has names ending in _L. Psi4 wants OFF" ON)
# CODE_COVERAGE should be tested and probably propagated to external projects
#option_with_flags(ENABLE_CODE_COVERAGE "Enables details on code coverage" OFF
#                  "-ftest-coverage")

######################### Process & Validate Options ###########################
include(autocmake_safeguards)
include(autocmake_omp)
include(autocmake_static_library)

################################# Main Project #################################
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(PN ${PROJECT_NAME})

# <<<  Build  >>>

set(raw_sources_list ai_disp.c balance.c bvec.c clapack.c disp.c
                     efp.c elec.c elec_terms.c int.c log.c
                     parse.c pol.c pol_direct.c stream.c swf.c
                     util.c xr.c)
set(src_prefix "src/")
string(REGEX REPLACE "([^;]+)" "${src_prefix}\\1" sources_list "${raw_sources_list}")

    # STATIC/SHARED on below governed by BUILD_SHARED_LIBS
add_library(efp ${sources_list})
set_target_properties(efp PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC}
                                     COMPILE_FLAGS "-std=c99")
set_source_files_properties(src/int.c PROPERTIES COMPILE_FLAGS -O1)
if(${BUILD_SHARED_LIBS})
    target_link_libraries(efp PRIVATE ${LIBC_INTERJECT})
endif()

if(NOT FRAGLIB_UNDERSCORE_L)
    file(GLOB _dotefps RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "fraglib/*.efp")
    foreach(_dotefp ${_dotefps})
        file(READ ${_dotefp} _strefpL)
        string(REPLACE "_L" "" _strefp ${_strefpL})
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${_dotefp} ${_strefp})
    endforeach()
endif()

# <<<  Install  >>>

if(FRAGLIB_UNDERSCORE_L)
    install(DIRECTORY fraglib
            DESTINATION ${CMAKE_INSTALL_DATADIR}/${PN})
else()
    install(FILES fraglib/makefp.inp
            DESTINATION ${CMAKE_INSTALL_DATADIR}/${PN})
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fraglib
            DESTINATION ${CMAKE_INSTALL_DATADIR}/${PN})
endif()

    # headers NOT namespace protected
install(FILES ${src_prefix}/efp.h 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS efp
        EXPORT "${PN}Targets"
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# <<<  Export Interface  >>>

target_compile_definitions(efp INTERFACE USING_${PN})
target_include_directories(efp INTERFACE
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# <<<  Export Config  >>>

set(CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/cmake/${PN}")
configure_package_config_file(${PN}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
                                 VERSION 1.3.0
                                 COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(EXPORT "${PN}Targets"
        NAMESPACE "${PN}::"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

# Notes on CMake: CMake builds libefp, the shared or static library,
#   but not efpmd because math detection would add a whole lot of extra code,
#   and Psi4 doesn't use it. Also not implemented is testing with CTest.

