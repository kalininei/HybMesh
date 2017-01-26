#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./hmmapping_test

set breakpoint pending on
b healgrid.cpp:251
run

source ../../src/libs/dbgscripts/skiplist
