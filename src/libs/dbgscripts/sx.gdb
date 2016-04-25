source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx fromdoc/intro_shark.py -silent
set breakpoint pending on
b contours.cpp:81
run

