@echo off
set HMCOM="C:\Program Files\HybMesh\bin\hybmesh.exe" -silent -sx

echo "============ doing unit tests"
%HMCOM% cont_test.py
%HMCOM% exclude_test.py
%HMCOM% bgrid_test.py
%HMCOM% unite_test.py
%HMCOM% imex_test.py
%HMCOM% mapping_test.py
%HMCOM% extrude_test.py

echo "============ doing doc tests"
%HMCOM% fromdoc/ex_setbtype.py
%HMCOM% fromdoc/ex_exclude.py
%HMCOM% fromdoc/ex_unite_grids.py
%HMCOM% fromdoc/ex_bgrid.py
%HMCOM% fromdoc/ex_partcontour.py
%HMCOM% fromdoc/ex_extrude.py
echo "--- intro_house"
%HMCOM% fromdoc/intro_house.py
echo "--- intro_pentagon"
%HMCOM% fromdoc/intro_pentagon.py
echo "--- intro_hcyl"
%HMCOM% fromdoc/intro_hcyl.py
echo "--- intro_shark"
%HMCOM% fromdoc/intro_shark.py
echo "--- intro_naca4415"
%HMCOM% fromdoc/intro_naca4415.py


pause