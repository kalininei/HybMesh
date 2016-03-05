.. _functionality:

Functionality
==============

Terminology
-----------

The program operates with two geometric object types: grid and contour.
Each contour and each grid has its own unique internal name which
is used for addressing.

Contour is a set of points connected in a defined order with each contour
section bearing a boundary type feature. Contour can be open or closed or
even multiply-connected (bounding a multiply connected domain).
Their direction could not be set implicitly. If fact
all contours which bound a domain apply left direction of traversal: 
all outer subcontours are traced counterclockwise and inner ones - 
in the opposite.

Grid is a set of cells which are defined as a sequences of points.
Grid cell can have arbitrary number of points,
the only restriction applied to a grid cell is that it can not
be multiply connected.

.. warning::
   
  Not every grid format supports arbitrary cells.
  Most FEM grids (like those used in GMsh) could contain
  only triangular or quadrangular cells. This should be
  taken into account during grid exporting.

Each grid contains its own bounding contour which is
referenced as a grid contour.
It includes all boundary nodes of the grid along with grid boundary features.
Most procedures which take contour as an invariable input parameter (e.g.
set boundary types, exclude contour from a grid etc.) could also be done 
using grid contours addressed by a grid internal name.

Boundary type is a non-geometric object which is defined
by integer positive boundary index (zero is the default value for non-defined
segments) and unique boundary name.
Grid exporting procedures try to use boundary names and indices
if export format supports it.


.. _gridimp:

Grid Imposition
---------------
This is the basic HybMesh operation. Generally it takes two independent
grids which have non-zero domain intersection and composes them into a single grid.
The domain of resulting grid is exactly equal to the domain of geometrical union of parent grids domains and
its cells reflect the original grids cells everywhere except for a zone around the line of parent grids
contact which is triangulated providing smooth cell size transition. This zone is later referenced as a
*buffer zone*.


Order of imposition matters. For clarity sake we call the first of two original grids the *base grid*
and the second one -- *imposed grid*. **Buffer is always built within the base grid**. Cells
of *imposed grid* are transfered to the resulting grid mostly untouched
(except for a few boundary vertices near grids intersection zone. See :ref:`fix-bnd-option` for details).
So, as you can see on the picture below, by swapping the grid roles we obtain different resulting grid geometry.

.. figure:: grid_imposition1.png
   :height: 600 px

   fig1. Basic imposition example

Hybmesh also supports imposition of the grids chain. In this regime a sequence of
imposition operations are performed over a list of grids. Each operation takes the result of previous operation
as a *base grid* and use the next grid in given list as an *imposed grid*.
You should carefully define the order of grids in a input list to get desirable result.


.. figure:: grid_imposition2.png
   :height: 400 px

   fig2. Chain imposition example


2. Imposition only takes place if parent grid has non-zero and non one-point intersection area.

3. Boundary features.

Buffer zone size
++++++++++++++++
Buffer zone size.

Sometimes it could be useful to compose grids without building triangulated buffer.

.. _fix-bnd-option:

Fix boundary nodes
++++++++++++++++++

Empty holes
+++++++++++



.. _bgrids:

Building Boundary Grids
-----------------------
TODO
