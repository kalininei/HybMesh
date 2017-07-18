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

   New python interface functions:

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

   New python interface functions:

   * :func:`extrude_grid`
   * :func:`partition_contour` 
   * :func:`export3d_grid_vtk`
   * :func:`export3d_grid_msh`
   * :func:`add_custom_rect_grid`

   Python interface functions with new options:

   * *periodic_pairs* for :func:`export_grid_msh`
* **0.4.1**:
   New algorithms:

   * build 3D grid by revolution of 2D grid around arbitrary `xy` vector
   * 2D/3D grid export to tecplot ascii \*.dat format

   New geometry constructors:
   
   * quadrangle grid in a circle area

   New python interface functions:

   * :func:`revolve_grid`
   * :func:`add_circ_rect_grid` 
   * :func:`export3d_grid_tecplot`
   * :func:`export_grid_tecplot`
   * :func:`export_contour_tecplot`

* **0.4.2**:
   New algorithms:

   * inverse Laplace algorithm for grid mapping,
   * Laplace and orthogonal algorithms for custom rectangular grid constructor,
   * Laplace and orthogonal algorithms for quadrangle grid in a circle area constructor.

   Python interface functions with new options:

   * *algo* = ``direct_laplace``, ``inverse_laplace`` for :func:`map_grid`,
   * *return_invalid* for :func:`map_grid`,
   * *algo* = ``laplace``, ``orthogonal_circ``, ``orthogonal_rect`` for :func:`add_circ_rect_grid`,
   * *algo* = ``direct_laplace``, ``inverse_laplace``, ``orthogonal`` for :func:`add_custom_rect_grid`.

* **0.4.3**:
   New algorithms:

   * linear and hermite transfinite interpolations
   * fixed number of edges for contour partition

   Python interface functions with new options:

   * *algo* = ``linear_tfi``, ``hermite_tfi`` for :func:`add_custom_rect_grid`,
   * *hermite_tfi_w* for :func:`add_custom_rect_grid`,
   * *nedges* for :func:`partition_contour`.

* **0.4.4**:
   New algorithms:

   * create contour using parametric splines

   New python interface functions:

   * :func:`create_spline_contour`

* **0.4.5**:
   New features:

   * unite grids with buffer filled with recombined mostly quadrangular mesh
   * constrained unstructured meshing with triangle, recombined and pebi grids
   * rectangle and radial grid prototypes with custom segmentation
   * procedure for concave cells removing
   * matched contour partition with 2D size conditions
   * numeric segment partition
   * native file formats introduced and described
   * export 3d grid to gmsh format
   * exporting list of 2d grids or contours to one file of any supported format
   * importing data from native files

   New geometry constructors:

   * stripe grid prototype
   * regular hexagonal grid prototype
   

   New python interface functions:

   * :func:`triangulate_domain`
   * :func:`pebi_fill` 
   * :func:`stripe`
   * :func:`add_unf_hex_grid`
   * :func:`partition_segment`
   * :func:`matched_partition`
   * :func:`export_all_hmd`
   * :func:`export3d_grid_gmsh`
   * :func:`export3d_grid_hmg`
   * :func:`import_grid_hmg`
   * :func:`import_contour_hmc`
   * :func:`import_all_hmd`
   * :func:`load_project`
   * :func:`registered_grids3d`


   Python interface functions with new options:

   * *buffer_fill* for :func:`unite_grids`
   * *custom_x, custom_y* for :func:`add_unf_rect_grid`
   * *custom_rads, custom_archs* for :func:`add_unf_circ_grid`
   * *crosses* for :func:`partition_contour`
   * *convex_cells* for :func:`heal_grid`

* **0.4.6**:
   New features:

   * forced connection of not connected open contours
     using shifting and stretching
   * decomposition of multiply connected edges into simple contours
   * subdivision of given contour into subcontours with user defined
     end points
   * building contour segmentation using recommended sizes defined by
     length or normalized length contour coordinates 
   * definition of target contour segment and fixed points for contour
     partition algorithm
   * mapping of grids into reflected domains
   * finding of an object point closest to the user defined site
   * calculation of 3D domain volume 
   * selection of contour closest to user defined point
   * extraction of grid surface
   * building tetrahedral mesh in closed 3D domain with arbitrary faces
   * exporting and importing of surfaces to native file format
   * merging of 3D grids with coincident boundaries


   New python interface functions:

   * :func:`grid3_bnd_to_surface`
   * :func:`decompose_contour`
   * :func:`extract_subcontours`
   * :func:`connect_subcontours`
   * :func:`merge_grids3`
   * :func:`info_grid3d`
   * :func:`info_surface`
   * :func:`get_point`
   * :func:`domain_volume`
   * :func:`pick_contour`
   * :func:`export3d_surface_hmc`
   * :func:`import3d_surface_hmc`
   * :func:`tetrahedral_fill`
   * :func:`registered_surfaces`


   Python interface functions with new options:

   * *keep_pts*, *start*, *end*, *algo* = ``ref_weights``, ``ref_lengths`` for
     :func:`partition_contour`
   * *is_reversed* for :func:`map_grid`
   * *return_invalid* for :func:`add_custom_rect_grid`


* **0.5.0**:
   New features:

   * Object oriented programming interfaces for C++, C#, Matlab, Java and Python languages
   * Windows installation as a Python2.7 package
   * Readonly access to geometric data of hybmesh objects
   * Grid to contour snapping algorithm
   * Consistent exceptions mechanism
   * More advanced standard output verbosity level control
   * Simplified function signatures for building boundary grids and grid union

   New python interface functions:

   * :func:`tab_cont2`, :func:`tab_grid2`, :func:`tab_surf3`, :func:`tab_grid3`
   * :func:`snap_grid_to_contour`
   * :func:`unite_grids1`, :func:`build_boundary_grid1`


* **0.5.1**:
   New python interface functions:

   * :func:`remove_cells`
   * :func:`inscribe_grid`
   * :func:`insert_grid_constraints`
