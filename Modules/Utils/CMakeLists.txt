project(Utils)

set(HEADER
    ssmUtils.h
    ssmOptionsBase.h
    ssmExtractionOptions.h
    ssmAlignmentOptions.h
	ssmReferenceOptions.h
)

add_custom_target(Utils ALL SOURCES ${HEADER})
