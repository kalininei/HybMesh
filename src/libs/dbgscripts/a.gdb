#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./hmblay_test

set breakpoint pending on
b unite_grids.cpp:231
run

# source ../../src/libs/dbgscripts/skiplist
