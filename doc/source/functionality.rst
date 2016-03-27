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
(except for a few boundary vertices near grids intersection zone. See :ref:`fixbnd` for details).
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


.. _fixbnd:

Fix boundary nodes
++++++++++++++++++

This option defines treatment of *base grid* boundary vertices if they get into buffer zone.
Sometimes building a smooth grid within the buffer demands remeshing those boundaries. However
this can lead to change of the initial domain area and loss of some boundary features. Therefore
user is given a choice whether to allow algorithm to move such vertices or not.

If *fix boundary nodes* option is set to True then all boundary vertices presented in *base grid* and
*overlaid grid* will be present in the result grid if they lie on its domain boundary.
With this option is on it is guaranteed that:

* *overlaid grid* is passed to result grid without any changes
* shape of domain intersection of input grids is precisely preserved
* boundary features of output grid exactly replicate input grid features

However there are some possible drawbacks of this option.
Picture below illustrates superposition of two square grid with complicated boundary set.
First example shows the result of operation with fixed boundary vertices.
Due to the fact that points of intersection don't hit any of the *overlaid grid*  
initial vertices two hanging nodes have appeared in the result.
Furthermore since some of boundary nodes of *base grid* lied too close to
these intersection points highly skewed triangles were built in the buffer zone.

Second example on the picture below shows the same operation without *fix boundary nodes* option.
In order to get smoother grid two vertices of *overlaid grid* were moved to intersection locations and
buffer zone boundary segments were remeshed. As a result we've completely lost blue and magenta boundary
segments but the resulting grid don't contain any hanging nodes or highly skewed cells.

.. figure:: grid_imposition7.png
   :width: 700 px

   fig6. Superposition with and without fixing boundary nodes


.. _zero-angle-app:
   
Zero angle approximation
++++++++++++++++++++++++

By default only boundary nodes which lie on a straight line (form an angle of 180 degree) could be moved when *fix boundary nodes* option
is off. This guaranties the exact preservation of domain intersection shape.
However if grids domain is formed by smooth curved lines the option *fix boundary nodes = False* will take no effect
since all points on such lines will be treated as corner points. Option *Zero angle approximation* (:math:`\alpha_0`) allows user
to define which boundary polyline turns should be considered negligible and be treated as straight angles.
With non-zero :math:`\alpha_0` all boundary vertices which lie within buffer zone and provide
turn between :math:`[180-\alpha_0, 180+\alpha_0]` will be considered as candidates for moving or removing.

The effect of :math:`\alpha_0` option is shown at figure 7. Both results here were obtained with *fix boundary nodes = False*.
The first was done with :math:`\alpha_0=0` hence all arc points were preserved and very coarse cell size transition occurred in the
bottom of the buffer zone. In the second example arc boundary segment of the buffer zone was remeshed to get better grid quality. However
due to loss of some shape forming nodes in the latter case result domain doesn't exactly equal
input grid intersection domain.

.. figure:: grid_imposition8.png
   :width: 700 px

   fig7. Superposition without fixing boundary nodes and different 
   zero angle :math:`{\alpha}_0` values

.. _emptyholes:

Empty holes
+++++++++++
If this option of grid superposition is set to true than all hulls presented at *overlaid grid* will be 
preserved as hulls in the result grid. Otherwise these hulls will be filled according to
general algorithm. The effect of this option is shown at picture below.

.. figure:: grid_imposition6.png
   :width: 600 px

   fig8. Superposition with/without *Empty holes* option


.. _bgrids:

Building Boundary Grids
-----------------------
TODO

.. _gridmappings:

Grid Mapping
------------
TODO
