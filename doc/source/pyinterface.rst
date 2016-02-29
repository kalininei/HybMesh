.. _pyinterf:

Python Interface
===================

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

List of Functions
-----------------
.. include:: functab

Object Procedures
------------------

These geometric procedures are used for both grid 
and contour objects addressed via string identifiers
returned by object constructors.

.. autofunction:: remove_geom

.. autofunction:: remove_all

.. autofunction:: move_geom

.. autofunction:: rotate_geom

.. autofunction:: scale_geom

.. autofunction:: copy_geom

Contour Procedures
------------------

Contour procedures could be done for both user and
grid contours. In the latter case grid identifier
instead of contour identifier should be used.

.. autofunction:: add_boundary_type

.. autofunction:: set_boundary_type

.. autofunction:: create_contour

.. autofunction:: add_rect_contour

.. autofunction:: add_circ_contour

.. autofunction:: grid_bnd_to_contour

.. autofunction:: simplify_contour

.. autofunction:: unite_contours

.. autofunction:: clip_domain

Grid Procedures
---------------

Prototypes
++++++++++

Building uniform grids in a primitive areas.

.. autofunction:: add_unf_rect_grid
.. autofunction:: add_unf_circ_grid
.. autofunction:: add_unf_ring_grid

Transformations
+++++++++++++++
See :ref:`functionality` for procedures details.

.. autofunction:: exclude_contours
.. autofunction:: unite_grids
.. autofunction:: build_boundary_grid

.. autoclass:: BoundaryGridOptions()
   :members: __init__, uniform_partition, incremental_partition

.. autofunction:: heal_grid

Information
-----------
Some general characteristics of program state and geometrical objects could be obtained by
set of special commands.

.. autofunction:: check_compatibility
.. autofunction:: info_grid
.. autofunction:: info_contour
.. autofunction:: registered_contours
.. autofunction:: registered_grids
.. autofunction:: registered_btypes
.. autofunction:: domain_area
.. autofunction:: skewness


Exports
-------
Possible Grid export formats are

* Vtk file format could be used to open a grid in `Paraview`.
  No boundary types are saved for this format.
  To see boundary types use separate vtk contour export command.
* :ref:`hmg-file`
* Fluent mesh format 
* GMsh file format

The current work flow can be also save to a :ref:`hmp-file` which
can be opened in a gui application.

.. autofunction:: export_grid_vtk
.. autofunction:: export_grid_hmg
.. autofunction:: export_grid_msh
.. autofunction:: export_grid_gmsh
.. autofunction:: export_contour_vtk
.. autofunction:: export_contour_hmc
.. autofunction:: save_project

Imports
-------

.. autofunction:: import_grid_msh
.. autofunction:: import_grid_gmsh
.. autofunction:: import_contour_ascii

Introductory Examples
---------------------
.. toctree::
   
  intro_example1
  intro_example2
  intro_example3
