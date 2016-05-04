#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./hmmaping_test
set breakpoint pending on
b rectangle_grid_builder.cpp:125
run 
