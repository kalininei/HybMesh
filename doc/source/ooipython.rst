Usage
^^^^^

Python wrapper contains a source file ``Hybmesh.py``
and shared library ``core_hmconnection_py``.
To use it in your application (either Python 2 or Python 3)
simply copy the source file into application directory and import
it as normal python file as

.. code-block:: py

  from Hybmesh import Hybmesh

Shared library can be left where it is.
Its absolute path along with absolute
path to installed hybmesh executable is
written into **Hybmesh.hybmesh_lib_path** and
**Hybmesh.hybmesh_exec_path** static properties respectively.
Those fields can be reassigned before first construction of
**Hybmesh** object if default values are not correct.

**Hybmesh** class supports creation under python ``with``
statement. By calling it as

.. code-block:: py

  from Hybmesh import Hybmesh

  with Hybmesh() as hm:
      # ... operations ...

it is guaranteed that superclass object will
be destroyed immediately at the end of ``with`` block and
hybmesh subprocess will be closed.
Otherwise the subprocess will wait until
garbage collector executes superclass destructor.


Introductory Example
^^^^^^^^^^^^^^^^^^^^
The following example illustrates numerical
computation of an integral of a function (Gaussian hill) on a grid.

It uses :ref:`circrect_grid` prototype grid with three different
outer step sizes. For each cell of the grid it calculates approximate
center point where it computes target function and multiplies the result by
the cell area.

Here paths to directories containing hybmesh executable and 
shared library are defined explicitly using relative paths.
If hybmesh was properly installed in the current system and Hybmesh.py file
is taken from the installation directory this could be omitted.

The following code could be run both on Python 2 and Python 3
interpreters.

.. literalinclude:: ../../testing/bindings/py/fromdoc/func_integral.py

