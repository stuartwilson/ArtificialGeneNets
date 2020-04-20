# AUTOMOC requires 2.8.6. Per-target include directories would require 2.8.8
cmake_minimum_required(VERSION 3.1)
project(ArtificialGeneNets)

message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  (This can be changed with `cmake -DCMAKE_INSTALL_PREFIX=/some/place`")

# From CMAKE_SYSTEM work out which of __OSX__, __GLN__, __NIX__ or
# __WIN__ are required
message(STATUS "Operating system: " ${CMAKE_SYSTEM})
if(CMAKE_SYSTEM MATCHES Linux.*)
set(MORPH_HOST_DEFINITION "-D__GLN__")
elseif(CMAKE_SYSTEM MATCHES BSD.*)
set(MORPH_HOST_DEFINITION "-D__NIX__")
elseif(APPLE)
set(MORPH_HOST_DEFINITION "-D__OSX__")
elseif(CMAKE_SYSTEM MATCHES Win.*)
set(MORPH_HOST_DEFINITION "-D__WIN__")
else()
message(ERROR "Operating system not supported: " ${CMAKE_SYSTEM})
endif()

set (CMAKE_CXX_STANDARD 17)

# Add the host definition to CXXFLAGS along with other switches
if (APPLE)
# OpenMP: According to
# https://stackoverflow.com/questions/39979836/using-openmp-with-c11-on-mac-os#39983933
# it should now be a case of:
#
# You should be able to get away without the libomp that's mentioned
# in that post, as I don't use any of the runtime library OpenMP
# features. However, the code did not compile without it, so we've moved the links to omp to the later apple if case.
#
set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -g -Xpreprocessor -fopenmp")
else()
# To use Intel compiler, call cmake as: cmake -DCMAKE_CXX_COMPILER=icpc ..
if (CMAKE_CXX_COMPILER_ID MATCHES Intel)
set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -g -xHOST -O3 -qopenmp -D__ICC__")
else() # GCC
set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -g -Wno-unused-result -Wno-unknown-pragmas -march=native -O3 -fopenmp")
endif()
endif()

# Be cavalier about comparing signed and unsigned ints ;)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")

# Important additional GL compiler flags.
#
# Use vendor implementations (not legacy softare impl.)
set(OpenGL_GL_PREFERENCE "GLVND")
# Important flags to set (perhaps these belong in the morphologica headers)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL3_PROTOTYPES -DGL_GLEXT_PROTOTYPES")

# Lib finding - we have to do this, so that the correct include paths
# are present for OpenGL, OpenCV, etc; they're all included by
# morphologica header files.
find_package(OpenCV REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(X11 REQUIRED)

if(APPLE)
find_package (HDF5 REQUIRED)
set (HDF5LIBLINK ${HDF5_C_LIBRARY_hdf5})
set (JSONLIBLINK /usr/local/lib/libjsoncpp.dylib)
include_directories(/opt/local/include/libomp)
link_libraries(-L/opt/local/lib/libomp omp)
else()
set (LIB_TYPE SHARED) # STATIC or SHARED
string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)
find_package (HDF5 NAMES hdf5 COMPONENTS C ${SEARCH_TYPE} REQUIRED)
set (HDF5LIBLINK ${HDF5_C_${LIB_TYPE}_LIBRARY})
set (JSONLIBLINK jsoncpp)
endif(APPLE)

find_package(glfw3 QUIET)

include_directories(${OpenCV_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR} ${X11_INC_SEARCH_PATH} ${HDF5_INCLUDE_DIR})

# Find libmorphologica using pkgconfig
find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
message(STATUS "Have pkg_config, searching for libmorphologica...")

if(APPLE)
# Manually set PKG_CONFIG_PATH on Mac, as mac ports with by default search only in /opt/local
set(ENV{PKG_CONFIG_PATH}  "ENV{PKG_CONFIG_PATH}:${CMAKE_INSTALL_PREFIX}/lib/pkgconfig")
endif(APPLE)

pkg_check_modules(MORPH REQUIRED libmorphologica)
if (MORPH_FOUND)
message(STATUS "pkg_config found libmorphologica")
include_directories(${MORPH_INCLUDEDIR})
include_directories(${MORPH_INCLUDEDIR}/morph)
else()
message(FATAL_ERROR "You need libmorphologica from github.com/ABRG_Models/morphologica")
endif(MORPH_FOUND)
else()
message(WARNING "There's no pkg-config on this system to check for libmorphologica. "
"You can set MORPH_LIB with `cmake -DMORPH_LIB=/pathto/lib` ...")
endif()

add_subdirectory(src)

# On Apple, have to explicitly add the X11 library directory
if(APPLE)
  link_directories(/usr/X11R6/lib)
endif()

# For debugging of variables:
set(DEBUG_VARIABLES OFF)
if(DEBUG_VARIABLES)
get_cmake_property(_variableNames VARIABLES)
foreach (_variableName ${_variableNames})
message(STATUS "${_variableName}=${${_variableName}}")
endforeach()
endif(DEBUG_VARIABLES)
