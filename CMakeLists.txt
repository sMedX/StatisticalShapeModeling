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

add_executable(maskToSurfaceConstructor MaskToSurfaceConstructor.cxx)
target_link_libraries(maskToSurfaceConstructor ${statismo_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(decimateSurface DecimateSurface.cxx)
target_link_libraries(decimateSurface ${statismo_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(referenceShapeConstructor ReferenceShapeConstructor.cxx)
target_link_libraries(referenceShapeConstructor ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(surfaceToImageRegistration SurfaceToImageRegistration.cxx)
target_link_libraries(surfaceToImageRegistration ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(shapeModelToSurfaceRegistration ShapeModelToSurfaceRegistration.cxx)
target_link_libraries(shapeModelToSurfaceRegistration ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(crossValidationTest CrossValidationTest.cxx)
target_link_libraries(crossValidationTest ${statismo_LIBRARIES} ${ITK_LIBRARIES})
