#build as archive since only hmmath is supposed to use it
add_library(${DSCPACK_TARGET} OBJECT dscpack.f )
set_target_properties(${DSCPACK_TARGET} PROPERTIES COMPILE_FLAGS "-cpp -fPIC")

