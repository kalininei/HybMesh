source ../../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../../src/py/hybmesh.py -sx intro_hcyl.py
set breakpoint pending on
b grid.cpp:425
# b buffergrid.cpp:232
run 

