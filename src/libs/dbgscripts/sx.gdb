source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx cont_test.py -silent
set breakpoint pending on
b hmcport.cpp:302
run 

