source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/HybMesh.py -sx bgrid_test.py
set breakpoint pending on
b bgrid_impose.cpp:297
run 

