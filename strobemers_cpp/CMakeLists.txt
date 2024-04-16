cmake_minimum_required(VERSION 3.22)

project(
        StrobeMap
        VERSION 0.1
        LANGUAGES CXX
)

#search openMP/ZLIB on system
find_package(OpenMP)
find_package(ZLIB)

# libraries present in the repo
add_library(
        kstream
        INTERFACE
            kseq++.hpp
            seqio.hpp
)
set_target_properties(kstream PROPERTIES LINKER_LANGUAGE CXX)

add_library(
        robinhood
        INTERFACE
            robin_hood.h
)
set_target_properties(robinhood PROPERTIES LINKER_LANGUAGE CXX)

################################
# set compilation options

set(GCC_COMPILE_OPTIONS "-Wall;-Wno-multichar;-Wno-comment;-Wno-sign-compare;-Wno-deprecated;-Wno-reorder;-Winline")
set(GCC_COMPILE_DEBUG_OPTIONS "${GCC_COMPILE_OPTIONS};-ggdb;-O0")
set(GCC_COMPILE_RELEASE_OPTIONS "${GCC_COMPILE_OPTIONS};-O3;-mavx2")

################################
# StrobeMap executable

add_executable(
        ${PROJECT_NAME}
        main.cpp
        index.cpp
)
# link
target_link_libraries(
        ${PROJECT_NAME}
        PUBLIC
            robinhood
            kstream
            OpenMP::OpenMP_CXX
            ZLIB::ZLIB
)
# set target dir for binary
set_target_properties(
        ${PROJECT_NAME}
        PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)
# set compilation options
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Debug>:${GCC_COMPILE_DEBUG_OPTIONS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Release>:${GCC_COMPILE_RELEASE_OPTIONS}>")

################################
# create a library, using index.hpp
# and make it searchable via find_package()

add_library(
        strobemer
        STATIC
)
target_sources(
        strobemer
        PRIVATE
            index.cpp
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/index.hpp>
            $<INSTALL_INTERFACE:include/index.hpp>
)

# the following is quite useless as both headers and sources ar in the same dir
# but kept for future refactoring and good cmake practices
target_include_directories(
        strobemer
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src
        PUBLIC
            # where top-level project will look for the library's public headers
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
            # where external projects will look for the library's public headers
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set_target_properties(
        strobemer
        PROPERTIES
        # debug postfix to differenciate debug and release builds
        DEBUG_POSTFIX "d"
        # set target dir for library
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
)
# set compilation options
target_compile_options(strobemer PUBLIC "$<$<CONFIG:Debug>:${GCC_COMPILE_DEBUG_OPTIONS}>")
target_compile_options(strobemer PUBLIC "$<$<CONFIG:Release>:${GCC_COMPILE_RELEASE_OPTIONS}>")

# by default install to /install
# note that it is not CMAKE_INSTALL_PREFIX we are checking here
if(DEFINED CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(
        STATUS
        "CMAKE_INSTALL_PREFIX is not set\n"
        "Default value: ${CMAKE_INSTALL_PREFIX}\n"
        "Will set it to ${CMAKE_SOURCE_DIR}/install"
    )
    set(CMAKE_INSTALL_PREFIX
        "${CMAKE_SOURCE_DIR}/install"
        CACHE PATH "Where the library will be installed to" FORCE
    )
else()
    message(
        STATUS
        "CMAKE_INSTALL_PREFIX was already set\n"
        "Current value: ${CMAKE_INSTALL_PREFIX}"
    )
endif()

# without it public headers won't get installed
set(public_headers index.hpp)
# note that ${public_headers} has to be in quotes
set_target_properties(strobemer PROPERTIES PUBLIC_HEADER "${public_headers}")

# debug postfix
set_target_properties(strobemer PROPERTIES DEBUG_POSTFIX "d")

# DESTINATION
# definitions of CMAKE_INSTALL_LIBDIR, CMAKE_INSTALL_INCLUDEDIR and others
include(GNUInstallDirs)
# install the target and create export-set
install(
    TARGETS strobemer
    EXPORT "strobemerTargets"
    # these get default values from GNUInstallDirs, no need to set them
    #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
    #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    # except for public headers, as we want them to be inside a library folder
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/strobemer # include/SomeLibrary
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
)

# generate and install export file
install(EXPORT "strobemerTargets"
    FILE "strobemerTargets.cmake"
    NAMESPACE ${namespace}::
    DESTINATION cmake
)

include(CMakePackageConfigHelpers)

# generate the version file for the config file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/strobemerConfigVersion.cmake"
    VERSION "${version}"
    COMPATIBILITY AnyNewerVersion
)

# create config file
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/strobemerConfig.cmake"
    INSTALL_DESTINATION cmake
)

# install config files
install(
    FILES
    "${CMAKE_CURRENT_BINARY_DIR}/strobemerConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/strobemerConfigVersion.cmake"
    DESTINATION cmake
)