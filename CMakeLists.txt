cmake_minimum_required(VERSION 3.0.0)
project(libMATA)

set(CMAKE_CXX_STANDARD 17)

# 3rd party modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# Enable fancy colours
set(CMAKE_COLOR_MAKEFILE ON)
# Enable verbose output for 'make'
# set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_MACOSX_RPATH 1)

##############################################################################
#                                     CCACHE
##############################################################################
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

##############################################################################
#                             OUTPUT CUSTOMIZATION
##############################################################################

# Build with Clang? -- currently does nothing, initialize with
#   $ CXX=clang++ make debug
# to enable Clang
option(USE_CLANG "build with clang" OFF)
# option(USE_CLANG "build with clang" ON)

##############################################################################
#                                DEPENDENCIES
##############################################################################

# find_package(Valgrind)

# Doxygen - voluntary
find_package(Doxygen)
if (DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
       ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
    )
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
else()
    message(WARNING "Doxygen not found. Documentation will not be built")
endif()

##############################################################################
#                               COMPILER FUN
##############################################################################

# enable C++
enable_language(CXX)

# get the git revision
include(GetGitRevisionDescription)
git_describe(GIT_DESCRIBE)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

message("-- Default C++ compiler: ${CMAKE_CXX_COMPILER}")

##############################################################################
#                                 PROJECT
##############################################################################

# Include CTest so that sophisticated testing can be done now
include(CTest)
enable_testing()

# if(CMAKE_COMPILER_IS_GNUCXX)
    include(CodeCoverage)
    setup_target_for_coverage(${PROJECT_NAME}_coverage tests coverage)
# endif()

##############################################################################
#                                SUB-MODULES
##############################################################################

# Process subdirectories
add_subdirectory(src)
add_subdirectory(cli)
add_subdirectory(3rdparty/re2)
add_subdirectory(3rdparty/simlib)
add_subdirectory(3rdparty/cudd)

install(TARGETS libmata ARCHIVE DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/mata DESTINATION include FILES_MATCHING PATTERN "*.hh")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/cudd/include/mata DESTINATION include FILES_MATCHING PATTERN "*.h*")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/simlib/include/mata DESTINATION include FILES_MATCHING PATTERN "*.hh")
