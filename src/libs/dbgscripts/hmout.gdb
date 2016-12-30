define contvtk
	if $argc ==1
		call HMCont2D::ECollection::SaveVtk($arg0, "_dbgout.vtk")
	else
		call HMCont2D::ECollection::SaveVtk($arg0, $arg1)
	end
end

define pntvtk
	if $argc ==1
		call HMCont2D::PCollection::SaveVtk($arg0, "_dbgout.vtk")
	else
		call HMCont2D::PCollection::SaveVtk($arg0, $arg1)
	end
end

define gridvtk
	if $argc ==1
		call GGeom::Debug::save_vtk($arg0, "_dbgout.vtk")
	else
		call GGeom::Debug::save_vtk($arg0, $arg1)
	end
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

define info_tree
	call HMCont2D::Debug::info_tree($arg0)
end

define info_ecollection
	call HMCont2D::Debug::info_ecollection($arg0)
end

define info_extpath
	call HMBlay::Debug::info_extpath($arg0)
end

define gg_contour
	call HMCont2D::Debug::geogebra_contour($arg0)
	if ($argc>1)
		call HMCont2D::Debug::geogebra_contour($arg1)
	end
	if ($argc>2)
		call HMCont2D::Debug::geogebra_contour($arg2)
	end
	if ($argc>3)
		call HMCont2D::Debug::geogebra_contour($arg3)
	end
	if ($argc>4)
		call HMCont2D::Debug::geogebra_contour($arg4)
	end
	if ($argc>5)
		call HMCont2D::Debug::geogebra_contour($arg5)
	end
	if ($argc>6)
		call HMCont2D::Debug::geogebra_contour($arg6)
	end
end

define gg_tree
	call HMCont2D::Debug::geogebra_tree($arg0)
	if ($argc>1)
		call HMCont2D::Debug::geogebra_tree($arg1)
	end
	if ($argc>2)
		call HMCont2D::Debug::geogebra_tree($arg2)
	end
	if ($argc>3)
		call HMCont2D::Debug::geogebra_tree($arg3)
	end
	if ($argc>4)
		call HMCont2D::Debug::geogebra_tree($arg4)
	end
	if ($argc>5)
		call HMCont2D::Debug::geogebra_tree($arg5)
	end
	if ($argc>6)
		call HMCont2D::Debug::geogebra_tree($arg6)
	end
end

define gg_ecol
	call HMCont2D::Debug::geogebra_ecollection($arg0)
	if ($argc>1)
		call HMCont2D::Debug::geogebra_ecollection($arg1)
	end
	if ($argc>2)
		call HMCont2D::Debug::geogebra_ecollection($arg2)
	end
	if ($argc>3)
		call HMCont2D::Debug::geogebra_ecollection($arg3)
	end
end


skip file /usr/include/c++/4.8/functional
skip file /usr/include/c++/4.8/bits/shared_ptr_base.h
skip file /usr/include/c++/4.8/bits/stl_algo.h
skip file /usr/include/c++/4.8/bits/stl_iterator.h
skip file /usr/include/c++/4.8/bits/allocator.h
skip file /usr/include/c++/4.8/bits/stl_vector.h
skip file /usr/include/c++/4.8/tuple
skip file /usr/include/c++/4.8/array
#skip file /usr/include/c++/4.8/bits/shared_ptr.h
