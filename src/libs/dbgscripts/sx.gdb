source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx extrude_test.py -silent
set breakpoint pending on
b construct_grid3d.cpp:694
b construct_grid3d.cpp:933
run

