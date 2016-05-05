cmake_minimum_required(VERSION 3.2)

project(StatisticalShapeModelling)

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

add_executable(binaryMask3DSurfaceConstructor BinaryMask3DSurfaceConstructor.cxx)
target_link_libraries(binaryMask3DSurfaceConstructor ${statismo_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(binaryMask3DProcessing BinaryMask3DProcessing.cxx)
target_link_libraries(binaryMask3DProcessing ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(shapeModelToImageRegistration ShapeModelToImageRegistration.cxx)
target_link_libraries(shapeModelToImageRegistration ${statismo_LIBRARIES} ${ITK_LIBRARIES})
