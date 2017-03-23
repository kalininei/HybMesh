Functions reference
-------------------

Object class
^^^^^^^^^^^^
This is a basic class for all geometrical objects.
It also provides methods for setting boundary features
and getting object structure information.

.. autoclass:: Hybmesh::Hybmesh.Object
    :members: free, dims, set_btypes_all, set_btypes, raw_vertices, raw_tab


Object2D class
^^^^^^^^^^^^^^
Basic class for 2D objects provides primitive geometric transformation routines.

.. autoclass:: Hybmesh::Hybmesh.Object2D
    :members: move, scale, rotate, reflect

Contour2D class
^^^^^^^^^^^^^^^
2D contour identifier.
See also :func:`Hybmesh.Hybmesh.Object2D`,
:func:`Hybmesh.Hybmesh.Object` for parent class members.

.. autoclass:: Hybmesh::Hybmesh.Contour2D
    :members: deepcopy, get_point, domain_area

Grid2D class
^^^^^^^^^^^^
2D grid identifier.
See also :func:`Hybmesh.Hybmesh.Object2D`,
:func:`Hybmesh.Hybmesh.Object` for parent class members.

.. autoclass:: Hybmesh::Hybmesh.Grid2D
    :members: deepcopy, skewness

Object3D class
^^^^^^^^^^^^^^
Basic class for 3D objects

.. autoclass:: Hybmesh::Hybmesh.Object3D
    :members: move, scale

Surface3D class
^^^^^^^^^^^^^^^
3D surface identifier.
See also :func:`Hybmesh.Hybmesh.Object3D`,
:func:`Hybmesh.Hybmesh.Object` for parent class members.

.. autoclass:: Hybmesh::Hybmesh.Surface3D
    :members: deepcopy, domain_volume


Grid3D class
^^^^^^^^^^^^
3D surface identifier.
See also :func:`Hybmesh.Hybmesh.Object3D`,
:func:`Hybmesh.Hybmesh.Object` for parent class members.

.. autoclass:: Hybmesh::Hybmesh.Grid3D
    :members: deepcopy

Hybmesh superclass
^^^^^^^^^^^^^^^^^^
.. rubric:: General operations

.. autoclass:: Hybmesh.Hybmesh
  :members: add_boundary_type,
     partition_segment,
     remove_all,
     pick_contour

.. rubric:: Callback manipulations

.. autoclass:: Hybmesh.Hybmesh
  :members:
    stdout_verbosity,
    assign_callback,
    reset_callback

.. rubric:: Contour operations

.. autoclass:: Hybmesh.Hybmesh
  :members:
    grid_bnd_to_contour,
    create_contour,
    create_spline_contour,
    extract_subcontours,
    add_rect_contour,
    add_circ_contour,
    simplify_contour,
    decompose_contour,
    unite_contours,
    clip_domain,
    connect_subcontours,
    partition_contour_const,
    partition_contour_ref_points,
    partition_contour_ref_lengths,
    matched_partition

.. rubric:: 2D grid operations

.. autoclass:: Hybmesh.Hybmesh
  :members:
    add_unf_rect_grid,
    add_unf_rect_grid1,
    add_unf_circ_grid,
    add_unf_ring_grid,
    add_unf_hex_grid_in_hex,
    add_unf_hex_grid_in_rect,
    add_triangle_grid,
    add_custom_rect_grid,
    add_custom_rect_grid_htfi,
    add_circ_rect_grid,
    stripe,
    triangulate_domain,
    pebi_fill,
    build_boundary_grid1,
    exclude_contours,
    unite_grids1,
    map_grid,
    heal_grid,
    snap_grid_to_contour


.. rubric:: Surface operations

.. autoclass:: Hybmesh.Hybmesh
  :members:
    grid3_bnd_to_surface

.. rubric:: 3D grid operations

.. autoclass:: Hybmesh.Hybmesh
  :members:
    tetrahedral_fill,
    extrude_grid,
    revolve_grid,
    merge_grids3

.. rubric:: Export

.. autoclass:: Hybmesh.Hybmesh
  :members:
    export_grid_vtk,
    export_grid_hmg,
    export_grid_msh,
    export_grid_gmsh,
    export_grid_tecplot,
    export3d_grid_vtk,
    export3d_grid_msh,
    export3d_grid_tecplot,
    export3d_grid_gmsh,
    export3d_grid_hmg,
    export_contour_vtk,
    export_contour_hmc,
    export_contour_tecplot,
    export3d_surface_hmc,
    export_all_hmd

.. rubric:: Import

.. autoclass:: Hybmesh.Hybmesh
  :members:
    import_grid_hmg,
    import_grid_msh,
    import_grid_gmsh,
    import3d_grid_hmg,
    import_contour_ascii,
    import_contour_hmc,
    import3d_surface_hmc

