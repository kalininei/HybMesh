#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./crossgrid_test
set breakpoint pending on
#b procgrid.cpp:433
b grid.cpp:530
run 
