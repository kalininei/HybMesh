define contvtk
	if $argc ==1
		call HMCont2D::ECollection::SaveVtk($arg0, "_dbgout.vtk")
	else
		call HMCont2D::ECollection::SaveVtk($arg0, $arg1)
	end
end

document contvtk
visualize HMCont2D data to _dbgout.vtk file or specified file
Usage: contvtk arg0 [arg1]
	arg0 -- HMCont2D::EdgeGeom object
	arg1 -- string filename
end

define ppvec
	if ($argc==1)
		set $i = 0
		set $sz = $arg0._M_impl._M_finish - $arg0._M_impl._M_start
	end
	if ($argc==2)
		set $i = $arg1
		set $sz = $arg0._M_impl._M_finish - $arg0._M_impl._M_start
	end
	if ($argc==3)
		set $i = $arg1
		set $sz = $arg2
	end
	while ($i<$sz)
		printf "[%i] = ", $i
		print **($arg0._M_impl._M_start + $i)
		set $i = $i + 1
	end
end

define info_contour
	call HMCont2D::Debug::info_contour($arg0)
end

define gg_contour
	call HMCont2D::Debug::geogebra_contour($arg0)
end


#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

file ./hmmath_test

skip file /usr/include/c++/4.8/functional
skip file /usr/include/c++/4.8/bits/shared_ptr_base.h
skip file /usr/include/c++/4.8/bits/stl_algo.h
skip file /usr/include/c++/4.8/bits/stl_iterator.h
skip file /usr/include/c++/4.8/bits/allocator.h
skip file /usr/include/c++/4.8/bits/stl_vector.h

b main
run 
b dscpack.f:1940

#====================================
#call from gui/HybMesh.py directory with HybMesh
# gdb -x ../libs/dbgscripts/a.gdb

#file python
#set args ./HybMesh.py ~/.HybMesh/_debug.hmp
#set breakpoint pending on
#b boundary_layer_grid.cpp:49
#run 
