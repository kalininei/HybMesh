source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/HybMesh.py -sx unite_test.py
set breakpoint pending on
b buffergrid.cpp:250
run 

