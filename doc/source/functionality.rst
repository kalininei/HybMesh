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
TODO


.. _bgrids:

Building Boundary Grids
-----------------------
TODO
