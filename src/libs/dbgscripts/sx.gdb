source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx bgrid_test.py -silent
set breakpoint pending on
b contour.cpp:390
#b bgrid_impose.cpp:52 if (i==49 && j==3 && k==9)
run 

