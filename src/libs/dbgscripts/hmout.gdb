define contvtk
	call HM2D::Debug::save_edges_vtk($arg0)
end

define vertvtk
	call HM2D::Debug::save_vertices_vtk($arg0)
end

define cellvtk
	call HM2D::Debug::save_cells_vtk($arg0)
end

define gridvtk
	call HM2D::Debug::save_grid_vtk($arg0)
end

define wfvtk
	call HM2D::Grid::Debug::save_wf_vtk($arg0)
end

define numerid
	call HM2D::Debug::numer_all($arg0)
end

define checkgrid
	if ($argc==1)
		call HM2D::Grid::Debug::report_grid_problems($arg0)
	end
	if ($argc==2)
		call HM2D::Grid::Debug::report_grid_problems($arg0, $arg1)
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
	call HM2D::Debug::info_contour($arg0)
end

define info_cell
	call HM2D::Debug::info_cell($arg0, $arg1)
end

define info_tree
	call HM2D::Debug::info_tree($arg0)
end

define info_extpath
	call HMBlay::Debug::info_extpath($arg0)
end


define gg_contour
	call HM2D::Debug::geogebra_contour($arg0)
	if ($argc>1)
		call HM2D::Debug::geogebra_contour($arg1)
	end
	if ($argc>2)
		call HM2D::Debug::geogebra_contour($arg2)
	end
	if ($argc>3)
		call HM2D::Debug::geogebra_contour($arg3)
	end
	if ($argc>4)
		call HM2D::Debug::geogebra_contour($arg4)
	end
	if ($argc>5)
		call HM2D::Debug::geogebra_contour($arg5)
	end
	if ($argc>6)
		call HM2D::Debug::geogebra_contour($arg6)
	end
end

define gg_tree
	call HM2D::Debug::geogebra_tree($arg0)
	if ($argc>1)
		call HM2D::Debug::geogebra_tree($arg1)
	end
	if ($argc>2)
		call HM2D::Debug::geogebra_tree($arg2)
	end
	if ($argc>3)
		call HM2D::Debug::geogebra_tree($arg3)
	end
	if ($argc>4)
		call HM2D::Debug::geogebra_tree($arg4)
	end
	if ($argc>5)
		call HM2D::Debug::geogebra_tree($arg5)
	end
	if ($argc>6)
		call HM2D::Debug::geogebra_tree($arg6)
	end
end

define gg_ecol
	call HM2D::Debug::geogebra_ecollection($arg0)
	if ($argc>1)
		call HM2D::Debug::geogebra_ecollection($arg1)
	end
	if ($argc>2)
		call HM2D::Debug::geogebra_ecollection($arg2)
	end
	if ($argc>3)
		call HM2D::Debug::geogebra_ecollection($arg3)
	end
end
