#!/usr/bin/env python
"Main executable"
import sys
import os.path


# version
class HybMeshVersion:
    def __init__(self, s):
        """initializes version for string of "1.2.3" type"""
        def __only_nums(n):
            return ''.join(k for k in n if k.isdigit())
        self.nums = map(__only_nums, s.split('.'))
        self.nums = map(int, [n for n in self.nums if len(n) > 0])

    def __str__(self):
        return '.'.join(map(str, self.nums))

    def str_compare(self, v2):
        '(str v2) -> [-1 or 0 or 1]. Compares self and v2 string'
        v2 = HybMeshVersion(v2)
        for i in range(min(len(self.nums), len(v2.nums))):
            if self.nums[i] < v2.nums[i]:
                return -1
            elif self.nums[i] > v2.nums[i]:
                return 1
        return 0
hybmesh_version = HybMeshVersion('0.0.3')
if len(sys.argv) == 2 and sys.argv[1] == '--version':
    print hybmesh_version
    sys.exit(0)

# config file: set homedir, libdir
if os.path.isfile('./config_hybmesh'):
    # tries current directory
    confpath = './config_hybmesh'
elif os.path.isfile('/etc/HybMesh/config_hybmesh'):
    # tries /etc for unix systems
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
gv.configure(confpath, hybmesh_version)

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
