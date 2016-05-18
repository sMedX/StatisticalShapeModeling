cmake_minimum_required(VERSION 3.2)

project(StatisticalShapeModeling)

# Find Statismo
find_package(statismo REQUIRED)

set(statismo_SOURCE_DIR CACHE PATH "")
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

add_executable(maskToSurfaceConstructor MaskToSurfaceConstructor.cxx)
target_link_libraries(maskToSurfaceConstructor ${statismo_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(surfaceToImageRegistration SurfaceToImageRegistration.cxx)
target_link_libraries(surfaceToImageRegistration ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(shapeModelToSurfaceRegistration ShapeModelToSurfaceRegistration.cxx)
target_link_libraries(shapeModelToSurfaceRegistration ${statismo_LIBRARIES} ${ITK_LIBRARIES})
