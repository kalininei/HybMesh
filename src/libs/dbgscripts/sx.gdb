source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx tmp.py -silent
set breakpoint pending on
b tetrahedral.cpp:62
run

