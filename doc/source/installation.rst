.. _installation:

Installation
============

Linux Installation
------------------

Dependencies
""""""""""""

The following packages (with their dependencies) should be already present in the system:

* gcc C++ and Fortran compilers. C++ compiler should support C++11 standard.
  Minimum tested version is 4.8.3.
* Python2.7  with `decorator package <https://pypi.python.org/pypi/decorator>`_  installed.
* CMake tool. Minimum tested version is 3.0.
* Blas, Lapack, SuiteSparse libraries with development headers.
* libxml2 with development headers.

Additionally you will need Octave and Java Development Kit to install
Octave and Java bindings respectively. Other wrappers
could be built without native tools. 

These are terminal commands which provide all necessary dependencies
on some systems:

.. rubric:: OpenSuse

Tested on clean newly installed `OpenSuse 42.2` with factory defaults.

.. code-block:: bash

   >>> sudo zypper in cmake gcc-c++ gcc-fortran blas-devel lapack-devel suitesparse-devel libxml2-devel

.. rubric:: Debian based

Tested on clean newly installed `Debian 8.0` and `Ubuntu 16.04` with ``sudo`` available.

.. code-block:: bash
   
   >>> sudo apt-get install g++ gfortran cmake libblas-dev liblapack-dev libsuitesparse-dev python-pip
   >>> sudo pip2 install decorator

You may also wish to install ``cmake-qt-gui`` to be able to run ``make edit_cache``
to edit configuration in GUI mode.


.. rubric:: Arch Linux

Tested on clean newly installed `Arch Linux 2017.03.01` with ``sudo`` available.

.. code-block:: bash
   
   >>> sudo pacman -S python2 python2-pip gcc-fortran blas lapack cmake suitesparse libxml2
   >>> sudo pip2 install decorator

Building
""""""""
Download latest HybMesh source archive `here 
<https://github.com/kalininei/hybmesh/releases/latest>`_, 
unpack it, open terminal in the source directory and execute

.. code-block:: bash

   >>> mkdir build && cd build
   >>> cmake .. -DCMAKE_BUILD_TYPE=Release
   >>> make -j8
   >>> sudo make install

After step 2 you may invoke ``make edit_cache`` to edit installation configuration if needed.
If your system already provides `Gmsh <gmsh.info>`_ shared library you
can force HybMesh to use it instead of embedded (somewhat obsolete) one
by assigning ``FIND_LIBGMSH`` variable.
To turn off certain wrapper installations use ``INCLUDE_***_BINDINGS`` flags.
Change ``CMAKE_INSTALL_PREFIX`` field to alter target installation path.

Run ``hybmesh -v`` to check if installation has been done properly.
This should return installed version of hybmesh.

If ``python-sphinx`` is available you can also build html documentation
by executing ``make`` from *HybMesh/doc* directory.

Summary
"""""""

As a result of the above procedure (with default ``CMAKE_INSTALL_PREFIX=/usr/local/``):

* ``hybmeshpack`` is installed as ordinary Python2 package,
* hybmesh executable is copied to ``/usr/local/hybmesh``,
* Wrappers are copied to ``/usr/local/include/hybmesh/*`` subdirectories,
* Additional libraries are installed into ``/usr/local/lib/hybmesh``.


Windows Installation
--------------------
To install HybMesh under Windows use installation distributive
which could be found `here <https://github.com/kalininei/hybmesh/releases/latest>`_.
Only 64bit platform is supported.

By default hybmesh core binaries, html documentation and set of
programming language wrappers will be installed into directory
defined by user. Main hybmesh executable *hybmesh.exe* will be
placed to *bin* subdirectory. Wrappers will be copied into *include* subdirectory. They already
include compiled dynamic libraries and are ready to use (of
course you will need respective runtimes). A shortcut to documentation main file *index.html*
will be placed to root installation folder.

User may also choose to install *hybmeshpack* as a regular Python package (that is turned off
by default).
This will allow to utilize hybmesh :ref:`script interface<pyinterf>` within python environment
omitting *hybmesh.exe* call.
Target system should provide a compatible 64bit Python2.7
with `decorator package <https://pypi.python.org/pypi/decorator>`_  installed.
That is not needed, however, if
you want to use python :ref:`programming interface<ooipython>` which also supports Python3
linking. Note that *hybmeshpack* package is not self-sufficient. It refers to
libraries located in *hybmesh/lib* directory with their absolute paths.
If Hybmesh program is deleted or moved this package would become non-operable.
Since package installation relies on python *distutils* setup procedure, uninstallation
of Hybmesh will not result in removing of *hybmeshpack*.
