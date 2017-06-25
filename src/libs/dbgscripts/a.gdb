#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./hmnumeric_test

set breakpoint pending on
b laplace_bem2d.cpp:102
run

# source ../../src/libs/dbgscripts/skiplist
