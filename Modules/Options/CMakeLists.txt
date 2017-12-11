project(Options)

set(HEADER
    ssmOptionsBase.h
    ssmExtractionOptions.h
    ssmAlignmentOptions.h
	ssmReferenceOptions.h
	ssmCorrespondenceOptions.h
	ssmModelBuildingOptions.h
)

add_custom_target(Options ALL SOURCES ${HEADER})
