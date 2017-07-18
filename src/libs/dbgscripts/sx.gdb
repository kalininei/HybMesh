source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx tmp.py -silent

#set args a.py
#set follow-fork-mode child

set breakpoint pending on
b sizefun.cpp:867
run

source ../../src/libs/dbgscripts/skiplist
