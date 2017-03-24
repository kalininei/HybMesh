.. HybMesh documentation master file, created by
   sphinx-quickstart on Mon Feb 15 08:43:06 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HybMesh Introduction
====================

HybMesh is a general-purpose grid generator
which implements locally-structured meshing approach.
With that approach grids are assembled from prebuilt structured grid segments
by superposition operations.
Grid segments are built using provided set of grid prototypes.
HybMesh includes algebraic mapping functionality which could be applied to
map those prototypes into any given non-regular area.

Generally, locally-structured grid generation workflow includes three basic steps:

* constructing grid prototypes,
* mapping prototypes into non-regular geometry (optionally),
* assembling the resulting grid by the superposition procedure.

Currently the above pattern is implemented for 2D grids generation only.
3D functionality includes building 3D grids from 2D grids using
extrusion and revolution procedures.

HybMesh could also be used to mesh arbitrary 2D and 3D domains
with fully unstructured grids. This functionality is provided
by internal calls of embedded `Gmsh <http://gmsh.info>`_ and
`TetGen <http://wias-berlin.de/software/tetgen/>`_ library routines respectively.

HybMesh application could be installed on **Windows** and **Linux** platforms.
Currently it provides only scripting interface based on Python2 language syntax.
Also a set of programming languages wrappers are included in installation
distributive, so HybMesh could be utilized as an embedded meshing engine
for **C++**, **C# (Mono)**, **Matlab (Octave)**, **Python 2/3**
and **Java** applications.

Example Grids
=============

+-----------------------------+-----------------------------+
|                                                           |
|.. figure:: overview_3.png                                 |
|   :height: 200 px                                         |
|   :figclass: align-center                                 |
|                                                           |
+-----------------------------+-----------------------------+
|                                                           |
|.. figure:: overview_2.png                                 |
|   :height: 200 px                                         |
|   :figclass: align-center                                 |
|                                                           |
+-----------------------------+-----------------------------+
|                             |                             |
|.. figure:: overview_1.png   | .. figure:: picintro_6_7.png|
|   :height: 300 px           |    :height: 300 px          |
|                             |                             |
+-----------------------------+-----------------------------+


Contents
========

.. toctree::
   :maxdepth: 2

   installation
   usage
   functionality
   prototypes
   pyinterface
   oointerface
   fileformats
   vers



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
