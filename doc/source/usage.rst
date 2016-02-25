
Program Usage
=============

.. _consoleapp:

Console Application
-------------------

Console application is a `hybmesh` or `hybmesh.exe` program which should
be invoked from the OS shell.
It is designed mostly for execution of hybmesh scripts (see :ref:`pyinterf`):

.. code-block:: console

  > hybmesh -sx yourscript.py [-silent]

Optional `-silent` flag is used to suppress all callback output from
the program. Possible stdout output from the user script will still be shown.

This program can also be used to execute project work flow from
the current till the last command (see :ref:`hmp-file`).

.. code-block:: console

  > hybmesh -s proj.hmp [-sgrid gname fmt filename] [-sproj filename.hmp] [-silent]

Here `proj.hmp` is an input work flow file.

`-sgrid gname fmt filename`: saves grid with internal name 'gname' to filename
using defined format (`hmg`: :ref:`hmg-file`, `vtk`: vtk format) at the end
of execution

`-sproj filename.hmp`: save work flow at the end of execution

`-silent`: suppress callbacks

Additional features of this program are

.. code-block:: console

  > hybmesh -help

prints usage help;

.. code-block:: console

  > hybmesh -v

prints installed hybmesh version;

.. code-block:: console

  > hybmesh -u

checks for latest release update on the github repository of the project
(Internet connection is required).


Gui
---
TODO

