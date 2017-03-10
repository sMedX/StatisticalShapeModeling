project(SSM)

set(HEADER
	ssmTypes.h
    ssmTransformInitializer.h
    ssmSurfaceToLevelSetImageFilter.h
    ssmSurfaceToLevelSetImageFilter.hxx
    ssmSurfaceToImageRegistrationMethod.h
    ssmSurfaceToImageRegistrationMethod.hxx
    ssmShapeModelToImageRegistrationMethod.h
    ssmShapeModelToImageRegistrationMethod.hxx
    ssmShapeModelToImageMetric.h
    ssmShapeModelToImageMetric.hxx
    ssmPointSetToPointSetMetrics.h
    ssmPointSetToImageMetrics.h
    ssmMeshPropertiesCalculator.h
    ssmMeshPropertiesCalculator.hxx
    ssmBinaryImageToLevelSetImageFilter.h
)

add_custom_target(SSM ALL SOURCES ${HEADER})
