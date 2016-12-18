source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx rrj2.py -silent
set breakpoint pending on
b grid.cpp:520 if buffer_size > 0
run

