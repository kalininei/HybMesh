#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./hmmaping_test
set breakpoint pending on
#b domapping.cpp:140
b femgrid43.cpp:100
run 
