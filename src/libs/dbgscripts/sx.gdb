source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
#set args ../../src/py/hybmesh.py -sx bugs/1/botscript1.py -silent
set args ../../src/py/hybmesh.py -sx tmp.py -silent
set breakpoint pending on
b trigrid.cpp:426
run

