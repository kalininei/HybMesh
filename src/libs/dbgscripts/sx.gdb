source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx fromdoc/intro_hcyl.py -silent
set breakpoint pending on
b bgrid.cpp:82
run

