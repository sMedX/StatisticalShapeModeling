project(SSM)

set(HEADER
    ssmTypes.h
	ssmImage3DMeshSource.h
	ssmImage3DMeshSource.hxx
    ssmInitializeSpatialTransform.h
	ssmShapeModelMultiTransform.h
	ssmShapeModelMultiTransform.hxx
    ssmMeshToLevelSetImageFilter.h
    ssmMeshToLevelSetImageFilter.hxx
    ssmMeshToImageRegistrationMethod.h
    ssmMeshToImageRegistrationMethod.hxx
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
