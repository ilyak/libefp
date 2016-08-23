## Libefp and CMake

#### Building

    cmake -H. -Bobjdir
    cd objdir && make
    make install

The build is also responsive to

- static/shared toggle `BUILD_SHARED_LIBS`
- install location `CMAKE_INSTALL_PREFIX`
- `name/name_L` in library fragments toggle `FRAGLIB_UNDERSCORE_L`
- `CMAKE_Fortran_COMPILER`, `CMAKE_C_COMPILER`, and `CMAKE_Fortran_FLAGS`

See [CMakeLists.txt](CMakeLists.txt) for options details and additional options.
All these build options should be passed as `cmake -DOPTION`.

#### Detecting

This project installs with `libefpConfig.cmake`, `libefpConfigVersion.cmake`, and `libefpTargets.cmake` files suitable for use with CMake [`find_package()`](https://cmake.org/cmake/help/v3.2/command/find_package.html) in `CONFIG` mode.

- `find_package(libefp)` - find any libefp libraries and headers
- `find_package(libefp 1.3.0 EXACT CONFIG REQUIRED COMPONENTS static)` - find libefp exactly version 1.3.0 built with static libraries; abort on failure

See [libefpConfig.cmake.in](libefpConfig.cmake.in) for details of how to detect the Config file and what CMake variables and targets are exported to your project.

#### Using in CMake-based project

After `find_package(libefp ...)`,

- test if package found with `if(${libefp_FOUND})` or `if(TARGET libefp::efp)`
- link to library (establishes dependency), including header and definitions configuration with `target_link_libraries(mytarget libefp::efp)`
- include header files using `target_include_directories(mytarget PRIVATE $<TARGET_PROPERTY:libefp::efp,INTERFACE_INCLUDE_DIRECTORIES>)`
- compile target applying `-DUSING_libefp` definition using `target_compile_definitions(mytarget PRIVATE $<TARGET_PROPERTY:libefp::efp,INTERFACE_COMPILE_DEFINITIONS>)`
