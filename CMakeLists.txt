cmake_minimum_required(VERSION 3.2)

project(StatisticalShapeModeling)

set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

# Find Boost
set(Boost_LIBRARY_DIR "" CACHE PATH "")
find_package(Boost 1.60 REQUIRED COMPONENTS filesystem program_options)
include_directories(${Boost_INCLUDE_DIRS})

# Boost setup
if(WIN32)
   # Disable autolinking in boost
   add_definitions(-DBOOST_ALL_NO_LIB)
   add_definitions(-DBOOST_ALL_DYN_LINK)
endif()

# Find Statismo
set(statismo_SOURCE_DIR CACHE PATH "")
find_package(statismo REQUIRED)

if(NOT statismo_SOURCE_DIR)
  message(FATAL_ERROR "'statismo_SOURCE_DIR' is required but not specified")
endif()

list(APPEND statismo_INCLUDE_DIRS ${statismo_SOURCE_DIR}/modules/ITK/cli)
include_directories(${statismo_INCLUDE_DIRS})

# Find ITK
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

if(${ITKVtkGlue_LOADED})
else()
  message(FATAL_ERROR "'ITKVtkGlue' module is required but it was not enabled in ITK '${ITK_DIR}'")
endif()

# Find VTK
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# ------------------------------------------------------------------------------
# Modules and applications
# ------------------------------------------------------------------------------
include_directories(Modules/SSM)
include_directories(Modules/Utils)

add_subdirectory(Modules/SSM)
add_subdirectory(Modules/Utils)
add_subdirectory(Applications)