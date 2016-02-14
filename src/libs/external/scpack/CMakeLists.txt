#build as archive since only hmmath is supposed to use it
add_library(${SCPACK_TARGET} OBJECT scpdbl.f sclibdbl.f)
set_target_properties(${SCPACK_TARGET} PROPERTIES COMPILE_FLAGS "-fPIC")
