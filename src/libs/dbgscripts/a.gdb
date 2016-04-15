#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

source ../../src/libs/dbgscripts/hmout.gdb

file ./hmgrid3d_test
set breakpoint pending on
b fluent_export_grid3d.cpp:62
run 
