C++ Bindings for HybMesh
========================

Usage
^^^^^

Unlike all other wrappers C++ interface contains only single header file ``Hybmesh.hpp``
without any other dependencies.
It could be copied into client application project directory
or included directed from its location.
To build that file C++ compiler should be set to support C++11 standard.

Since there is no dynamically linked libraries only a
directory containing Hybmesh executable should be defined.
Default path is set in wrapper source file through

.. code-block:: cpp

  #define DEFAULT_HYBMESH_EXE_PATH "/hybmesh/install/directory/"

macro directive.

To use another path create Hybmesh instance with
optional string argument:

.. code-block:: cpp

  Hybmesh hm("custom/hm/path/");

Besides basic geometrical and exception classes
Hybmesh superclass also provides 2 additional nested classes defining 2D
and 3D Point which are used to pass point arguments to hybmesh methods.

All geometric nested classes (including point classes) provide
special static method ``None()`` which returns a respective
object that is treated as ``None-object``.
It is used for those hybmesh methods which accept such arguments
(for example **right** and **top** contour arguments
of :func:`Hybmesh.Hybmesh.add_custom_rect_grid`)

For detailed description of all methods see
:ref:`python wrapper reference<ooifuncref>`
and documentation in ``Hybmesh.hpp`` file.

.. _cppintro:

Helloworld Example
^^^^^^^^^^^^^^^^^^

After hybmesh installation create a directory with copied
*Hybmesh.hpp* and the following *test.cpp* file.

.. code-block:: cpp
    :caption: test.cpp

    #include <iostream>
    #include "Hybmesh.hpp"
    
    int main(){
    	Hybmesh hm;
    	Hybmesh::Grid2D g2 = hm.add_unf_rect_grid(
    		Hybmesh::Point2(0, 0), Hybmesh::Point2(1, 1), 2, 2);
    	std::cout<<"number of cells: "<<g2.dims()[2]<<std::endl;
    }

With gcc (or MinGW) compiler open terminal at created directory and
run the following terminal commands 

.. code-block:: bash

    >>> g++ -o test test.cpp -std=c++11
    >>> ./test

With VisualStudio compiler run Developer Command Prompt and invoke

.. code-block:: bash

    >>> cl /EHsc test.cpp
    >>> test.exe


Introductory Example
^^^^^^^^^^^^^^^^^^^^

In the following example a unite operation between
regular rectangular grid and ring grid is performed.
Two grids are built by different Hybmesh instances
in two different threads in parallel.

The secondary thread builds ring grid and saves it to file.
The first thread builds rectangular grid, waits until
secondary thread finishes its job, reads ring grid from file
and performs grid union.

.. literalinclude:: ../code_preproc_out-intro_multithreaded.cpp
  :language: c++
  :tab-width: 4
