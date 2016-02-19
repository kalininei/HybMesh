#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./hybmesh_contours2d_test
b main
b hybmesh_contours2d_test.cpp:260
run 
