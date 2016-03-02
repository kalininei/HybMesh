source ../../src/libs/dbgscripts/hmout.gdb

#call from testing directory
file python
set args ../../src/py/hybmesh.py -sx fromdoc/intro_hcyl.py -silent
set breakpoint pending on
b algos.cpp:172
#b canonic_bgrid.cpp:93 if ( ((double(*)(double))fabs)(h - 0.0113038749)<1e-5  && bottom.lenght()>0.055 && bottom.length()<0.056 )
b canonic_bgrid.cpp:93 if (bottom.length()<0.056 )
b canonic_bgrid.cpp:203
run 

