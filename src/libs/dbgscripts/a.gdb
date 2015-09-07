define contvtk
	if $argc ==1
		call HMCont2D::SaveVtk($arg0, "_dbgout.vtk")
	else
		call HMCont2D::SaveVtk($arg0, $arg1)
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

document ppvec
prints pointer vector.
Usage ppvec arg0 [start_index] [endindex]
end


#call from build/bin directory with crossgrid_test:
# gdb -x ../../src/libs/dbgscripts/a.gdb

file ./hybmesh_contours2d_test
#file ./crossgrid_test
#file ./hmblay_test
b main
run 
b hybmesh_contours2d_test.cpp:169

#====================================
#call from gui/HybMesh.py directory with HybMesh
# gdb -x ../libs/dbgscripts/a.gdb

#file python
#set args ./HybMesh.py ~/.HybMesh/_debug.hmp
#set breakpoint pending on
#b boundary_layer_grid.cpp:49
#run 
