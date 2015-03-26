#!/usr/bin/env python
" Main executable "

import sys
from PyQt4 import QtGui
from PyQt4.QtCore import QT_VERSION_STR
from vtk import vtkVersion

# print used libraries versions
print "Python version:", ".".join(map(str, sys.version_info[0:3]))
print "Qt version:", QT_VERSION_STR
print "VTK version:", vtkVersion.GetVTKVersion()

# -- initialize qt application
app = QtGui.QApplication(sys.argv)

import globvars

# modes
if globvars.prog_options.debug_save_before:
    print "+++++ Debug save to: %s" % globvars.prog_options.debug_save_fn

# show main window
globvars.mainWindow.show()

# Vtk initialization after main window show().
# Otherwise vtk fails.
globvars.mainWindow.vtkWidget.iren().Initialize()

# open file from argument string
if len(sys.argv) > 1:
    globvars.Flows.xml_load(sys.argv[1])
else:
    globvars.actual_data().view_update()

# start gui loop
sys.exit(app.exec_())
