#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./hmgrid3d_test
set breakpoint pending on
b export_grid3d.cpp:165
#b construct_grid3d.cpp:122
run 
