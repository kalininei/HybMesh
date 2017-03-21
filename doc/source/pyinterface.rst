.. _pyinterf:

Python Script Interface
=======================

Script Invocation
-----------------
HybMesh python script is a normal `python2.7` file which
supports all capabilities of python language and can be invoked using both
HybMesh standalone application and python interpreter.
The latter demands HybMesh be installed as a python module
(see :ref:`installation`).

To invoke a script from python interpreter use a command

.. code-block:: console

   > python yourscript.py

For standalone :ref:`consoleapp` in Windows type

.. code-block:: bat

   > c:\path\to\bin\hybmesh.exe -sx yourscript.py

In Linux the default application install path `/usr/local/bin`
should present in your system paths, so there is no need to type full path
to executable.

.. code-block:: console

   > hybmesh -sx yourscript.py

HybMesh functionality is scattered across different modules of
:mod:`hybmeshpack` package. However all imports essential for scripting
are gathered in ``hybmeshpack.hmscripts.__init__`` module.
So the HybMesh script file should start with an import line

.. code-block:: python

   from hybmeshpack import hmscript

.. py:module:: hybmeshpack.hmscript

Design
------

All hybmesh operable geometrical objects
are represented by their internal unique string
identifiers. To get the full list of all registered identifiers
use ``'registered_*'`` functions.
If an identifier leaves current scope
it doesn't lead to disposing of respective object.
You should use :func:`remove_geom` function to explicitly
remove an object. Use :func:`remove_all_but` to
clean the garbage and leave only those objects 
which are actually needed.

To get full information about object geometry
(vertices coordinates and index based connectivity tables) use
``'tab_*'`` set of functions.
For performance reasons these functions return ctypes plain 1D arrays.
Use ordinary ``[i]`` notation to access entries of these arrays
or convert it to python list by ``pylist = ctypes_array[:]`` line.

Boundary (zone) types are defined by non-negative integers.
Default zone type (equals ``0``) is used if it was not explicitly defined.
You could use :func:`add_boundary_type`
function to link this integer to a string which will be
used in object export procedures.

Exception Handling
------------------

All hybmesh function may throw two types of exceptions:

.. autoexception:: ExecError
.. autoexception:: InvalidArgument

Before execution of any hybmesh function a quick
input data check is performed. It examines
argument types, ranges and consistency.
If it fails than :func:`InvalidArgument` is raised.
If output verbosity level is not set to zero then
detailed information regarding this error including
unity based index of invalid argument, its actual value and program expectations
will be flushed to console.

If program fails during function execution
then :func:`ExecError` will be raised.
Information about function arguments and traceback
will be printed to console.

Both exceptions could be safely caught and handled.


List of Functions
-----------------
.. include:: functab

General Procedures
------------------

.. autofunction:: check_compatibility
.. autofunction:: registered_contours
.. autofunction:: registered_grids
.. autofunction:: registered_surfaces
.. autofunction:: registered_grids3d
.. autofunction:: registered_btypes
.. autofunction:: remove_all
.. autofunction:: remove_all_but
.. autofunction:: copy_geom
.. autofunction:: move_geom
.. autofunction:: rotate_geom
.. autofunction:: scale_geom
.. autofunction:: reflect_geom
.. autofunction:: remove_geom
.. autofunction:: add_boundary_type
.. autofunction:: set_boundary_type
.. autofunction:: partition_segment


Contour Operations
------------------

.. autofunction:: info_contour
.. autofunction:: tab_cont2
.. autofunction:: domain_area
.. autofunction:: get_point
.. autofunction:: pick_contour
.. autofunction:: create_contour
.. autofunction:: create_spline_contour
.. autofunction:: add_rect_contour
.. autofunction:: add_circ_contour
.. autofunction:: add_circ_contour2
.. autofunction:: add_circ_contour3
.. autofunction:: unite_contours
.. autofunction:: clip_domain
.. autofunction:: simplify_contour
.. autofunction:: grid_bnd_to_contour
.. autofunction:: partition_contour
.. autofunction:: matched_partition
.. autofunction:: decompose_contour
.. autofunction:: extract_subcontours
.. autofunction:: connect_subcontours

2D Grid Operations
------------------

.. autofunction:: info_grid
.. autofunction:: tab_grid2
.. autofunction:: skewness
.. autofunction:: heal_grid
.. autofunction:: add_unf_rect_grid
.. autofunction:: add_unf_circ_grid
.. autofunction:: add_unf_ring_grid
.. autofunction:: add_unf_hex_grid
.. autofunction:: add_triangle_grid
.. autofunction:: add_custom_rect_grid
.. autofunction:: add_circ_rect_grid
.. autofunction:: stripe
.. autofunction:: triangulate_domain
.. autofunction:: pebi_fill
.. autofunction:: unite_grids
.. autofunction:: unite_grids1
.. autoclass:: BoundaryGridOptions()
   :members: __init__, uniform_partition, incremental_partition
.. autofunction:: build_boundary_grid
.. autofunction:: build_boundary_grid1
.. autofunction:: exclude_contours
.. autofunction:: map_grid


Surface Operations
------------------

.. autofunction:: info_surface
.. autofunction:: tab_surf3
.. autofunction:: grid3_bnd_to_surface
.. autofunction:: domain_volume

3D Grid Operations
------------------

.. autofunction:: info_grid3d
.. autofunction:: tab_grid3
.. autofunction:: extrude_grid
.. autofunction:: revolve_grid
.. autofunction:: tetrahedral_fill
.. autofunction:: merge_grids3

Export
------

.. autofunction:: export_grid_hmg
.. autofunction:: export_grid_vtk
.. autofunction:: export_grid_msh
.. autofunction:: export_grid_gmsh
.. autofunction:: export_grid_tecplot
.. autofunction:: export_contour_hmc
.. autofunction:: export_contour_vtk
.. autofunction:: export_contour_tecplot
.. autofunction:: export3d_grid_hmg
.. autofunction:: export3d_grid_vtk
.. autofunction:: export3d_grid_gmsh
.. autofunction:: export3d_grid_msh
.. autofunction:: export3d_grid_tecplot
.. autofunction:: export3d_surface_hmc
.. autofunction:: export_all_hmd
.. autofunction:: save_project

Import
------

.. autofunction:: import_grid_hmg
.. autofunction:: import_grid_msh
.. autofunction:: import_grid_gmsh
.. autofunction:: import_contour_hmc
.. autofunction:: import3d_grid_hmg
.. autofunction:: import3d_surface_hmc
.. autofunction:: import_all_hmd
.. autofunction:: load_project


Introductory Examples
---------------------
.. toctree::

  intro_example1
  intro_example2
  intro_example3
  intro_example4
  intro_example5
  intro_example6

.. toctree::
  :hidden:

  examples_include_files
