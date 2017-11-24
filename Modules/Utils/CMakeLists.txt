project(Utils)

set(HEADER
    ssmUtils.h
    ssmBaseOptions.h
    ssmSurfaceExtractionOptions.h
    ssmAlignmentOptions.h
	ssmReferenceOptions.h
)

add_custom_target(Utils ALL SOURCES ${HEADER})
