source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx fromdoc/intro_shark.py -silent
set breakpoint pending on
#b femgrid43.cpp:105
b  grid.cpp:467
run

