project(SSM)

set(HEADER
    ssmTypes.h
	ssmBinaryMask3DMeshSource.h
	ssmBinaryMask3DMeshSource.hxx
    ssmInitializeSpatialTransform.h
	ssmShapeModelMultiTransform.h
	ssmShapeModelMultiTransform.hxx
    ssmSurfaceToLevelSetImageFilter.h
    ssmSurfaceToLevelSetImageFilter.hxx
    ssmSurfaceToImageRegistrationMethod.h
    ssmSurfaceToImageRegistrationMethod.hxx
    ssmShapeModelRegistrationMethodBase.h
    ssmShapeModelRegistrationMethodBase.hxx
    ssmShapeModelToImageRegistrationMethod.h
    ssmShapeModelToImageRegistrationMethod.hxx
    ssmShapeModelToLevelSetImageMetric.h
    ssmShapeModelToLevelSetImageMetric.hxx
    ssmPointSetToPointSetMetrics.h
    ssmPointSetToImageMetrics.h
    ssmMeshPropertiesCalculator.h
    ssmMeshPropertiesCalculator.hxx
    ssmBinaryImageToLevelSetImageFilter.h
)

add_custom_target(SSM ALL SOURCES ${HEADER})
