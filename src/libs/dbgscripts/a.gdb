#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./hmblay_test
set breakpoint pending on
b main
run 
b hmblay_test.cpp:48
b grid.cpp: 430
