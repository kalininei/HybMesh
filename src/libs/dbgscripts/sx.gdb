source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx test.py -silent

#set breakpoint pending on
#b rectangle_grid_builder.cpp:262
#run

source ../../src/libs/dbgscripts/skiplist
