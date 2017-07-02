source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx others/rrj2.py -silent

#set args a.py
#set follow-fork-mode child

set breakpoint pending on
b cport_cont2d.cpp:674
run

source ../../src/libs/dbgscripts/skiplist
