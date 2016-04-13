source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx mapping_test.py -silent
set breakpoint pending on
#b  cport_cont2d.cpp:89
b rectangle_grid_builder.cpp:218
run

