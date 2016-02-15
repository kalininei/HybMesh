
Program Usage
=============

.. _consoleapp:

Console Application
-------------------

Console application is a `HybMesh`/`HybMesh.exe` program which should
be invoked from the os shell.
It is designed mostly for execution of HybMesh scripts (see :ref:`pyinterf`):

.. code-block:: console

  > HybMesh -sx yourscript.py [-silent]

Optional `-silent` flag is used to supress all callback output from
the program. Possible stdout output from the user script will still be shown.

This program can also be used to execute project workflow from
the current till the last command (see :ref:`hmp-file`).

.. code-block:: console

  > HybMesh -s proj.hmp [-sgrid gname fmt filename] [-sproj filename.hmp] [-silent]

Here `proj.hmp` is an input workflow file.

`-sgrid gname fmt filename`: saves grid with internal name gname to filename
using defined format (`hmg`: :ref:`hmg-file`, `vtk`: vtk format) at the end
of execution

`-sproj filename.hmp`: save workflow at the end of execution

`-silent`: supress callbacks

Additional features of this program are

.. code-block:: console

  > HybMesh -help

prints usage help;

.. code-block:: console

  > HybMesh -v

prints installed HybMesh version;

.. code-block:: console

  > HybMesh -u

checks for latest release update on the github repository of the project
(Internet connection is required).


Gui
---
TODO

