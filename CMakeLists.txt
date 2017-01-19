cmake_minimum_required(VERSION 3.2)

project(StatisticalShapeModeling)

# Find Statismo
find_package(statismo REQUIRED)

set(statismo_SOURCE_DIR CACHE PATH "")

if(NOT statismo_SOURCE_DIR)
  message(FATAL_ERROR "'statismo_SOURCE_DIR' is required but not specified")
endif()

list(APPEND statismo_INCLUDE_DIRS ${statismo_SOURCE_DIR}/modules/ITK/cli)

include_directories(${statismo_INCLUDE_DIRS})

# Boost setup
if(WIN32)
  # Disable autolinking in boost
  add_definitions(-DBOOST_ALL_NO_LIB)
endif()

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

add_executable(imageToSurface ImageToSurface.cxx)
target_link_libraries(imageToSurface ${statismo_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(alignmentOfSurfaces AlignmentOfSurfaces.cxx)
target_link_libraries(alignmentOfSurfaces ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(ssm-establish-correspondences ssm-establish-correspondences.cxx)
target_link_libraries(ssm-establish-correspondences ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(computeMeanSurface ComputeMeanSurface.cxx)
target_link_libraries(computeMeanSurface ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(shapeModelToSurfaceRegistration ShapeModelToSurfaceRegistration.cxx)
target_link_libraries(shapeModelToSurfaceRegistration ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(crossValidationTest CrossValidationTest.cxx)
target_link_libraries(crossValidationTest ${statismo_LIBRARIES} ${ITK_LIBRARIES})
