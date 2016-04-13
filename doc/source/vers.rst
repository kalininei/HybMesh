Versions History
=================

.. py:module:: hybmeshpack.hmscript

* **0.2.1**:
   basic python script interface introduced
* **0.3.0**:
   New algorithms:

   * grid mapping
   * zero angle approximation for grid superposition
   * geometry reflection over a straight line
   * vertex and corner projection types were added for boundary
     source contour start/end point definition

   New geometry constructors:
   
   * grid in a triangle area
   * circle contour defined by three points and
     by two points with curvature coefficient

   New python interface function:

   * :func:`remove_all_but`
   * :func:`reflect_geom` 
   * :func:`add_circ_contour2`
   * :func:`add_circ_contour3`
   * :func:`add_triangle_grid`
   * :func:`map_grid`

   Python interface functions with new options:

   * *zero_angle_approx* for :func:`unite_grids`
   * *project_to* for :class:`BoundaryGridOptions`
* **0.4.0**:
   New algorithms:

   * extrude 2D grid in z direction
   * 3D grid export to vtk and fluent msh formats
   * partition contour with constant and reference point defined step
   * fluent export with given periodic conditions
   * build rectangular grid on the basis of two or four open contours

   New python interface function:

   * :func:`extrude_grid`
   * :func:`partition_contour` 
   * :func:`export3d_grid_vtk`
   * :func:`export3d_grid_msh`
   * :func:`add_custom_rect_grid`

   Python interface functions with new options:

   * *periodic_pairs* for :func:`export_grid_msh`
