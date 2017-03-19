.. _oointerfaces:

High-Level Programming Interfaces
=================================

To use programming interface hybmesh should
be properly installed.
Needed files are located in ``include/*`` subdirectories of install directory.
By default on Linux systems it is ``/usr/local/include/hybmesh/``,
on Windows systems - ``C:\Program Files\Hybmesh\include\``.

The set of interface files for each language includes
a single native source file and (for all except C++) a C
shared library file. Source file is cross-platform
whereas library file is system dependent.

While being activated in target application,
those interfaces launch
the hybmesh executable as a background subprocess and
establish pipe communication between target (client) application
and hybmesh (server). Each interface command is translated into
hybmesh low-level python command and sent to
server. Then server response is parsed and 
translated back into native language data type.

Necessary paths to the hybmesh executable and respective
shared library are hardcoded into interface source files
during installation process. If for some reason
those default values are not correct (for example 
if the program was transfered from another system) they
could be changed either by manual edit of those files
or at runtime. See respective language chapters
and introductory examples for details.

Source files include very brief documentation on
their public methods. It contains method summary
and a link to the low-level interface function which
is used. To get detailed information on any needed function
find respective entry at :ref:`pyinterf` section.

General Design
--------------

Hybmesh high-level object orienting programming interface provides
a superclass dubbed **Hybmesh**.
While being created, an object of this class establishes
server connection.
On destruction it breaks the connection and all
derived objects become invalid.
Objects of this class should not be deep copied (cloned).

**Hybmesh** class contains a set of nested classes representing hybmesh
geometrical data objects which are nothing but
short string identifiers. Hence these objects can be cheaply
copied/passed by value, safely casted etc.
Couple of (mostly informational) methods
are implemented as methods of these classes however
main functionality is provided by the superclass methods.

The UML class diagram representing hybmesh interface is presented
below.

.. figure:: figs/oointerfaces_uml.png
   :width: 1200 px

All data representing object geometry is stored on server side
but could be accessed through :func:`Hybmesh.Hybmesh.raw_vertices`
and :func:`Hybmesh.Hybmesh.raw_tab` object methods.

Note that if an object is destructed (disposed) in terms
of target language it does not lead to object
destruction at server side. To remove object from hybmesh
operating list use :func:`Hybmesh.Hybmesh.Object.free` member function.
However even *free()* call does not guarantee release of resources.
The most reliable way to free memory resources is to
build all necessary grids, copy them to native side structures
and destroy the whole superclass object.


Multiple processes
^^^^^^^^^^^^^^^^^^

It is fine to have multiple instances of **Hybmesh** class at a time.
But the only way geometrical objects derived 
from one superclass can communicate with another
superclass is through file export/import mechanism.
Introductory :ref:`example<cppintro>` for C++ language illustrates
how two hybmesh subprocesses can be invoked in
parallel to build a single grid.

Exception Handling
^^^^^^^^^^^^^^^^^^

All interfaces except for Matlab support
exception mechanism. Two exception
classes are declared within **Hybmesh** superclass:
**ERuntimeError** and **EUserInterrupt**.
The first one is thrown when any error occurred
in hybmesh function. Detailed information about an
error could be extracted using exception methods.
The second exception is thrown if cancellation
was requested from callback function.
Matlab wrapper simply calls **error(msg)** function
in both cases.

Callbacks
^^^^^^^^^

All interfaces include callback mechanism with process progress reports
and cancellation support for time consuming operations.
Callback is a client side function with following general signature
(details depend on wrapper language)

.. code-block:: cpp

  int callback(string n1, string n2, double p1, double p2);

**n1** and **n2** are string names of current process and subprocess;
**p1** and **p2** are process and subprocess progress given in ``[0, 1]`` range.
**p2** could equal ``-1`` which means that subprocess reports are not supported
for current operation. Callback function normally returns ``0``;
for cancellation request it should return ``1``, in that case **EUserInterrupt**
exception will be raised and should be properly handled at the client side.

To assign certain callback function to hybmesh instance use
:func:`Hybmesh.Hybmesh.assign_callback` method. To reset it to default silent mode
use :func:`Hybmesh.Hybmesh.reset_callback`.

:ref:`Java<javaintro>` and :ref:`C#<csintro>` examples
use hybmesh to build simple
GUI applications which use :ref:`map grid<gridmappings>` and
:ref:`unite_grid<gridimp>`
hybmesh functionality respectively supplemented
with graphic callback.

By default standard console output of hybmesh is turned off
but it could be activated using
:func:`Hybmesh.Hybmesh.stdout_verbosity` method.

Functions reference
-------------------

All public method signatures present in source files 
are the same regardless the language. In fact
they were generated automatically from
the same pattern. Therefore in present documentation
only Python wrapper methods are described.

.. toctree::
   
   ooifuncref


Python
------

.. toctree::
   
   ooipython

C++
---

.. toctree::
   
   ooicpp

Java
----

.. toctree::
   
   ooijava

Matlab (Octave)
---------------

.. toctree::
   
   ooimatlab

C# (Mono)
---------

.. toctree::
   
   ooics
