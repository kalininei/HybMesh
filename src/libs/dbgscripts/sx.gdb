source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx fromdoc/intro_hcyl.py -silent

#set args a.py
#set follow-fork-mode child

set breakpoint pending on
b unite_grids.cpp:165
run

source ../../src/libs/dbgscripts/skiplist
