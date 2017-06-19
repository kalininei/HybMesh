source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx others/botscript1.py -silent

#set args a.py
#set follow-fork-mode child

set breakpoint pending on
b unite_grids.cpp:197
run

source ../../src/libs/dbgscripts/skiplist
