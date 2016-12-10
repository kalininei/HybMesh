source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx rrj2.py -silent
set breakpoint pending on
b cont_assembler.cpp:414
run

