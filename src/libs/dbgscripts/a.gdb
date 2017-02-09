#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./hmblay_test

set breakpoint pending on
b canonic_bgrid.cpp:526
run

source ../../src/libs/dbgscripts/skiplist
