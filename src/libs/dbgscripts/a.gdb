#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./hmblay_test

set breakpoint pending on
b infogrid.cpp:519 if nnn==8
#b hmblay_test.cpp:619
run



source ../../src/libs/dbgscripts/skiplist
