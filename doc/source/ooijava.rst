Java Bindings for HybMesh
=========================

Usage
^^^^^

Java wrapper consists of ``Hybmesh.java`` source file
and ``core_hmconnection_java`` library.
To compile it you need Java7 (or higher) development kit.
Source file can be copied to target application source directory.
There is no mechanism to define location of shared library
within java code. To make application see that library
it should be run with ``-Djava.library.path`` flag as

.. code-block:: bat

  java -Djava.library.path=directory/containing/core_hmconnection_java Application

To set custom path of hybmesh executable assign static property

.. code-block:: java

   Hybmesh.hybmesh_exec_path = "directory/containing/hybmesh";

before the first use of **Hybmesh** class.

Besides basic geometrical and exception classes
Hybmesh superclass also provides 2 additional nested classes defining 2D
and 3D Point which are used to pass point arguments to hybmesh methods.

Unfortunately Java lacks default function argument notation
which is heavily used by Hybmesh design.
So all arguments should be defined making function calls somewhat bulky.

Default argument values which are normally used
by other interfaces are listed in respective javadocs.
For non-primitive data types (i.e. geometrical objects,
arrays, strings etc.) arguments with default values can be assigned
with ``null``. If so they will be automatically reassigned
to their default values.

For example this is a function for importing surface as it is
defined in ``Hybmesh.java`` file.

.. code-block:: java

   /**
    * Imports surface from hybmesh native format.
    * 
    * See details in hybmeshpack.hmscript.import3d_surface_hmc().
    *
    * @param surfname if null then ""
    * @param allsurfs default is false
    */
   public Surface3D[] import3DSurfaceHmc(String fname, String surfname, boolean allsurfs)
           throws Hybmesh.EUserInterrupt, Hybmesh.ERuntimeError{
       // ....
   }

Here the **fname** argument should be defined by user, **surfname** could be
assigned with ``null`` (which will lead to default empty string),
**allsurfs** should be explicitly set to ``false`` to provide default behavior.

**Hybmesh** class implements *Autocloseable* interface so
superclass instances can be created within
try-with-resources block to guarantee
subprocess exit at its end. Otherwise it won't be closed
until it is captured by garbage collector.

.. code-block:: java

  try (Hybmesh hm = new Hybmesh()){
      //....
  }

For detailed description of all methods consult
:ref:`python wrapper reference<ooifuncref>`
and embedded javadocs of ``Hybmesh.java`` file.

Helloworld Example
^^^^^^^^^^^^^^^^^^

After installation of hybmesh program copy *Hybmesh.java* and the following
*Test.java* into empty directory.

.. code-block:: java
     :caption: Test.java

     class Test{
         public static void main(String[] args) throws Exception{
             try (Hybmesh hm = new Hybmesh()){
                 Hybmesh.Grid2D g2 = hm.addUnfRectGrid(
                         new Hybmesh.Point2(0, 0),
                         new Hybmesh.Point2(1, 1),
                         2, 2, null);
                 System.out.printf("number of cells: %d\n", g2.dims()[2]);
             }
         }
     }

Open terminal at this directory and execute (changing directory path to correct one)

.. code-block:: bash

    >>> javac Test.java Hybmesh.java
    >>> java -Djava.library.path=/directory/containing/core_hmconnection_java Test

.. _javaintro:

Introductory Example
^^^^^^^^^^^^^^^^^^^^

This example illustrates the creation of toy GUI
application which makes :ref:`mapping<gridmappings>` of regular unit square
grid into user-defined quadrangle with additional custom mapping points.
Graphics is written using Java swing library.

To compile and run this application create a directory
containing ``App.java`` and ``GridData.java`` files listed below along
with wrapper ``Hybmesh.java`` source file.
Then execute following commands (adjusting directory path)

.. code-block:: bash

    javac App.java GridData.java Hybmesh.java
    java -Djava.library.path=/path/to/core_hmconnection_java App

Here is the screenshot of the program as it appears in KDE system:

.. figure:: figs/javaexample.png
   :width: 600 px

To define custom mapping points turn on *Add refpoints* button
and click on the contours.
By turning on *Move refpoints* you can drag defined points
along contours and move target contour base vertices with the mouse.
With *Remove refpoints* click on defined points to delete them.
Basic grid segmentation is defined in a dialog activated by *Basic grid* button.
Grid mapping is performed at *Run mapping* click.

All defined points are indexed. The number and order of points
in basic and target contours should match as this is a requirement
of Grid mapping algorithm.

In this application Grid mapping procedure provides creation
of Progress graphic dialog with cancellation support.

Program source consists of two files:

* App.java contains visual interface code,
* GridData.java contains code needed to build and store hybmesh grids.

Procedure call with graphical callback support is provided
by *App.java.ProgressBarExecutor.execute()* static method.

.. Warning::

  For illustration and testing purposes callback function of current
  application provides intentional 0.3 second delay which
  is defined in *App.java:ProgressBarExecutor.Popup.callback()*
  function. It could be safely removed.


.. literalinclude:: ../code_preproc_out-App.java
    :language: java
    :caption: App.java
    :tab-width: 4

.. literalinclude:: ../code_preproc_out-GridData.java
    :language: java
    :caption: GridData.java
    :tab-width: 4


