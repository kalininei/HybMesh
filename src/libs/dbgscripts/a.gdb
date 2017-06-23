#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./crossgrid_test

set breakpoint pending on
b inscribe_grid.cpp: 20
run

# source ../../src/libs/dbgscripts/skiplist
