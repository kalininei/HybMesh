.. _example1:

Example 1
=========

This example illustrates usage of HybMesh algorithm for meshing
domain shown in :ref:`introfig1`

.. _introfig1:

.. figure:: picintro1.png
   :height: 400 px

   fig1. Meshing area

The building strategy is:

* build a set of primitive rectangular grids
* exclude pentagon (:ref:`introfig2`)
* unite all grids (:ref:`introfig3`)
* build a boundary grid and unite it with the previous result (:ref:`introfig4`)

Here is the script:

.. literalinclude:: ../../testing/py/fromdoc/intro_pentagon.py
   :end-before: ^^^^^^^^^^^^^^^

+--------------------------+--------------------------+--------------------------+
|.. _introfig2:            | .. _introfig3:           | .. _introfig4:           |
|                          |                          |                          |
|.. figure:: picintro2.png | .. figure:: picintro3.png| .. figure:: picintro4.png|
|   :height: 400 px        |    :height: 400 px       |    :height: 400 px       |
|                          |                          |                          |
|   fig2. Primitive grids  |    fig3. Superposition   |    fig4. Final result    |
+--------------------------+--------------------------+--------------------------+
