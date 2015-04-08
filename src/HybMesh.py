#!/usr/bin/env python
" Main executable "
import sys
import os.path

# version
hybmesh_version = '0.0.2'
if len(sys.argv) == 2 and sys.argv[1] == '--version':
    print hybmesh_version
    sys.exit(0)

# config file: set homedir, libdir
if os.path.isfile('./config_hybmesh'):
    confpath = './config_hybmesh'
elif os.path.isfile('/etc/HybMesh/config_hybmesh'):
    confpath = '/etc/HybMesh/config_hybmesh'
else:
    raise Exception('Config file not found')


from PyQt4 import QtGui
from PyQt4.QtCore import QT_VERSION_STR
from vtk import vtkVersion

# print used libraries versions
print "Python version:", ".".join(map(str, sys.version_info[0:3]))
print "Qt version:", QT_VERSION_STR
print "VTK version:", vtkVersion.GetVTKVersion()

# -- initialize qt application
app = QtGui.QApplication(sys.argv)

import HybMesh_gui.globvars as gv

# read configure
gv.configure(confpath)
if gv.prog_options.debug_save_before:
    print ">>> Debug save to %s" % gv.prog_options.debug_save_fn

# show main window
gv.mainWindow.show()

# Vtk initialization after main window show().
# Otherwise vtk fails.
gv.mainWindow.vtkWidget.iren().Initialize()

# open file from argument string
if len(sys.argv) > 1:
    gv.Flows.xml_load(sys.argv[1])
else:
    gv.actual_data().view_update()

# start gui loop
sys.exit(app.exec_())
