source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx test/test.py -silent
set breakpoint pending on
b grid.cpp:489
run

