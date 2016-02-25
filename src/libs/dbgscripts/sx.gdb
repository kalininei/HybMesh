source ../../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../../src/py/hybmesh.py -sx intro_house.py
set breakpoint pending on
b grid.cpp:493
run 

