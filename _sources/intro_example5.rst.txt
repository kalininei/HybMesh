.. _example5:

Example 5
=========

| 3D grid around an airfoil with circular vortex generator is
  built by extrusion of 2D grid in z direction.
| Airfoil data is loaded from external file (:ref:`naca4415_airfoil`).
| Resulting grid is saved to fluent \*.msh format with
  periodical merging of top-bottom and right-left surfaces.

.. figure:: picintro11.png
   :width: 700 px

Resulting 3D grid

------------------------------------

.. figure:: picintro9.png
   :width: 700 px

2D grid

-------------------------------------

.. figure:: picintro10.png
   :width: 700 px

2D grid in detail

-------------------------------------

.. literalinclude:: ../../testing/py/fromdoc/intro_naca4415.py
   :end-before: ^^^^^^^^^^^^^^^
