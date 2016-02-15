.. _pyinterf:

Python Interface
===================

Script Invocation
-----------------
HybMesh python script is a normal python 2.7 file which
supports all capabilities of python language and can be invoked using both
HybMesh standalone application and python interpreter.
The latter demands HybMesh be installed as a python module
(see :ref:`installation`).

To invoke a script from python interpreter use a command

.. code-block:: console

   > python yourscript.py

For standalone :ref:`consoleapp` in Windows system type

.. code-block:: bat

   > c:\path\to\bin\HybMesh.exe -sx yourscript.py

In Linux the default application install path `/usr/local/bin`
should present in your system paths, so there is no need to type full path
to executable.

.. code-block:: console

   > HybMesh -sx yourscript.py

HybMesh functionality is scattered across different modules of
:mod:`HybMeshPyPack` package. However all imports essential for scripting
are gathered in :mod:`HybMeshPyPack.hmscripts.__init__` module.
So the HybMesh script file should start with an import line

.. code-block:: python

   from HybMeshPyPack import hmscript

.. py:currentmodule:: HybMeshPyPack.hmscript

List of Functions
-----------------
* :func:`RemoveGeom`
* :func:`MoveGeom`
* :func:`RotateGeom`
* :func:`ScaleGeom`
* :func:`CopyGeom`
* :func:`AddBoundaryType`
* :func:`CreateContour`
* :func:`AddRectCont`
* :func:`SetBTypeToContour`
* :func:`GridBndToCont`
* :func:`SimplifyContour`
* :func:`UniteContours`
* :func:`AddUnfRectGrid`
* :func:`AddUnfCircGrid`
* :func:`AddUnfRingGrid`
* :func:`ExcludeContours`
* :func:`UniteGrids`
* :func:`BuildBoundaryGrid`
* :class:`BoundaryGridOption`
* :func:`ExportVTK`
* :func:`ExportHMG`
* :func:`ExportMSH`
* :func:`ExportGMSH`
* :func:`SaveProject`
* :func:`ImportGridMSH`
* :func:`ImportGridGMSH`
* :func:`ImportContourASCII`

Object Procedures
------------------

These geometric procedures are used for both grid 
and contour objects addressed via string identifiers
returned by object constructors.

.. autofunction:: HybMeshPyPack.hmscript.RemoveGeom

.. autofunction:: HybMeshPyPack.hmscript.MoveGeom

.. autofunction:: HybMeshPyPack.hmscript.RotateGeom

.. autofunction:: HybMeshPyPack.hmscript.ScaleGeom

.. autofunction:: HybMeshPyPack.hmscript.CopyGeom

Contour Procedures
------------------

Contour procedures could be done for both user and
grid contours. In the latter case grid identifier
instead of contour identifier should be used.

.. autofunction:: HybMeshPyPack.hmscript.AddBoundaryType

.. autofunction:: HybMeshPyPack.hmscript.CreateContour

.. autofunction:: HybMeshPyPack.hmscript.AddRectCont

.. autofunction:: HybMeshPyPack.hmscript.SetBTypeToContour

.. autofunction:: HybMeshPyPack.hmscript.GridBndToCont

.. autofunction:: HybMeshPyPack.hmscript.SimplifyContour

.. autofunction:: HybMeshPyPack.hmscript.UniteContours

Grid Procedures
---------------

Prototypes
++++++++++

Building uniform grids in a primitive areas.

.. autofunction:: HybMeshPyPack.hmscript.AddUnfRectGrid
.. autofunction:: HybMeshPyPack.hmscript.AddUnfCircGrid
.. autofunction:: HybMeshPyPack.hmscript.AddUnfRingGrid

Transformations
+++++++++++++++
See :ref:`functionality` for procedures details.

.. autofunction:: HybMeshPyPack.hmscript.ExcludeContours
.. autofunction:: HybMeshPyPack.hmscript.UniteGrids
.. autofunction:: HybMeshPyPack.hmscript.BuildBoundaryGrid


.. autoclass:: HybMeshPyPack.hmscript.BoundaryGridOption()
   :members: __init__

Exports
-------
Possible Grid export formats are

* Vtk file format could be used to open a grid in `Paraview`.
  No boundary types are saved for this format.
* :ref:`hmg-file`
* Fluent mesh format 
* GMsh file format

Current workflow can be also save to a :ref:`hmp-file` which
can be opened in a gui application.

.. autofunction:: HybMeshPyPack.hmscript.ExportVTK
.. autofunction:: HybMeshPyPack.hmscript.ExportHMG
.. autofunction:: HybMeshPyPack.hmscript.ExportMSH
.. autofunction:: HybMeshPyPack.hmscript.ExportGMSH
.. autofunction:: HybMeshPyPack.hmscript.SaveProject

Imports
-------

.. autofunction:: HybMeshPyPack.hmscript.ImportGridMSH
.. autofunction:: HybMeshPyPack.hmscript.ImportGridGMSH
.. autofunction:: HybMeshPyPack.hmscript.ImportContourASCII

Introductory Example
--------------------
This example illustrates HybMesh algorithm for meshing
domain shown in :ref:`introfig1`

.. _introfig1:

.. figure:: picintro1.png
   :height: 400 px

   fig1. Area for meshing

The building strategy is:

* build a set of primitive rectangular grids
* exclude pentagon (:ref:`introfig2`)
* unite all grids (:ref:`introfig3`)
* to get rid of oblong cells nearby the pentagon build a boundary grid and
  unite it with the result of previous imposition (:ref:`introfig4`)

Here is the script:

.. literalinclude:: ../../examples/intro.py

+--------------------------+--------------------------+--------------------------+
|.. _introfig2:            | .. _introfig3:           | .. _introfig4:           |
|                          |                          |                          |
|.. figure:: picintro2.png | .. figure:: picintro3.png| .. figure:: picintro4.png|
|   :height: 400 px        |    :height: 400 px       |    :height: 400 px       |
|                          |                          |                          |
|   fig2. Primitive grids  |    fig3. After imposition|    fig4. Final result    |
+--------------------------+--------------------------+--------------------------+
