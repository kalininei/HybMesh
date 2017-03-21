.. _installation:

Installation
============

The program can be installed either as a python package or as a standalone
application.

Linux Installation
------------------
TODO

Windows Installation
--------------------
To install HybMesh under Windows use installation distributive
which could be found at https://github.com/kalininei/hybmesh/releases/latest.
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
Target system should provide a compatible 64bit Python2 (>2.7)
with `decorator package <https://pypi.python.org/pypi/decorator>`_  installed.
That is not needed, however, if
you want to use python :ref:`programming interface<ooipython>` which also supports Python3
linking. Note that *hybmeshpack* package is not self-sufficient. It refers to
libraries located in *hybmesh/lib* directory with their absolute paths.
If Hybmesh program is deleted or moved this package would become non-operable.
Since package installation relies on python *distutils* setup procedure, uninstallation
of Hybmesh will not result in removing of *hybmeshpack*.
