project(SSM)

add_executable(ssm-build-shapemodel ssm-build-shapemodel.cxx)
target_link_libraries(ssm-build-shapemodel ${Boost_PROGRAM_OPTIONS_LIBRARY} ${statismo_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(ssm-extract-surface ssm-extract-surface.cxx)
target_link_libraries(ssm-extract-surface ${Boost_PROGRAM_OPTIONS_LIBRARY} ${statismo_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(ssm-align-surfaces ssm-align-surfaces.cxx)
target_link_libraries(ssm-align-surfaces ${Boost_PROGRAM_OPTIONS_LIBRARY} ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(ssm-build-reference ssm-build-reference.cxx)
target_link_libraries(ssm-build-reference ${Boost_PROGRAM_OPTIONS_LIBRARY} ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(ssm-establish-correspondence ssm-establish-correspondence.cxx)
target_link_libraries(ssm-establish-correspondence ${Boost_PROGRAM_OPTIONS_LIBRARY} ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(ssm-cross-validation ssm-cross-validation.cxx)
target_link_libraries(ssm-cross-validation ${Boost_PROGRAM_OPTIONS_LIBRARY} ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(ssm-fit-shapemodel-to-surface ssm-fit-shapemodel-to-surface.cxx)
target_link_libraries(ssm-fit-shapemodel-to-surface ${Boost_PROGRAM_OPTIONS_LIBRARY} ${statismo_LIBRARIES} ${ITK_LIBRARIES})
add_executable(ssm-specificity ssm-specificity.cxx)
target_link_libraries(ssm-specificity ${Boost_PROGRAM_OPTIONS_LIBRARY} ${statismo_LIBRARIES} ${ITK_LIBRARIES})

add_executable(ssm-fit-shapemodel-to-image ssm-fit-shapemodel-to-image.cxx)
target_link_libraries(ssm-fit-shapemodel-to-image ${Boost_PROGRAM_OPTIONS_LIBRARY} ${statismo_LIBRARIES} ${ITK_LIBRARIES})
