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
Their direction could not be set implicitly. In fact
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

Grid Superposition
------------------
This is the basic HybMesh operation. Generally it takes two independent
grids which have non-zero domain intersection and composes them into a single grid.
The domain of resulting grid is exactly equal to the domain of geometrical union of parent grid domains and
its cells reflect the original grids cells everywhere except for a zone around the line of parent grids
contact which is triangulated providing smooth cell size transition. This zone is later referenced as a
*buffer zone*.


Order of superposition matters. For clarity sake we call the first of two original grids the *base grid*
and the second one -- *overlaid grid*. **Buffer is always built within the base grid**. Cells
of *overlaid grid* are transfered to the resulting grid mostly untouched
(except for a few boundary vertices near grids intersection zone. See :ref:`fix-bnd-option` for details).
So, as you can see on the picture below, by swapping the grid roles we obtain different resulting grid geometry.

.. figure:: grid_imposition1.png
   :height: 600 px

   fig1. Basic superposition example

Hybmesh also supports superposition of grid chain. In this regime a sequence of
superposition operations are performed over a list of grids. Each operation takes the result of previous one
as a *base grid* and use the next grid in given list as an *overlaid grid*.
You should carefully define the order of grids in a input list to get desirable result.
On a picture below you can see the superposition result depending on given grids order.

.. figure:: grid_imposition2.png
   :height: 800 px

   fig2. Chain superposition example


Superposition with building a buffer grid only takes place if parent grid has non-zero and non single point intersection area.
Different operation results depending on relative position of input grids are presented below.

.. figure:: grid_imposition3.png
   :width: 600 px

   fig3. Superposition depending on types of given grid intersections.

If grid domains have no proper intersections (two last examples on the picture above) then the
resulting grid will contain cells from both given grids assembled to a single connectivity table.

Boundary features of superposed grid contour reflect boundary features of given grids.
If any boundary segment is contained in both *base grid* and *overlaid grid* then priority
will be given to features from the latter.


Buffer zone size
++++++++++++++++
Buffer zone is constructed as an area of all *base grid* cells which contain a vertex located
no further than given buffer zone size from contact line. Larger buffer zone provides
smoother triangle grid within the buffer (see picture below).

.. figure:: grid_imposition4.png
   :width: 600 px

   fig4. Superposition with different buffer sizes

Sometimes it is useful to superpose grids without building triangulated buffer. This could be
done by setting zero buffer zone size.
However, if vertices of *base grid* and *overlaid grid* do not coincide at contact line superposed
grid will contain hanging nodes (see second example at picture below). 
The necessity of superposing grids with zero buffer could be arisen e.g. while connecting grid to
an outer boundary grid built from its own contour.

.. figure:: grid_imposition5.png
   :width: 600 px

   fig5. Superposition with zero buffer size


.. _fix-bnd-option:

Fix boundary nodes
++++++++++++++++++
.. figure:: grid_imposition7.png
   :width: 700 px

   fig6. Superposition with/without *Empty holes* option

Zero angle approximation
++++++++++++++++++++++++
TODO

Empty holes
+++++++++++
If this option of grid superposition is set to true than all hulls presented at *overlaid grid* will be 
preserved as hulls in the result grid. Otherwise these hulls will be filled according to
general algorithm. The effect of this option is shown at picture below.

.. figure:: grid_imposition6.png
   :width: 600 px

   fig6. Superposition with/without *Empty holes* option



.. _bgrids:

Building Boundary Grids
-----------------------
TODO

.. _gridmappings:

Grid Mapping
------------
TODO
