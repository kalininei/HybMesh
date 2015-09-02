#!/usr/bin/env python
from PyQt4 import QtGui
from PyQt4.QtCore import QT_VERSION_STR
from vtk import vtkVersion

def main():
    # print used libraries versions
    print "Python version:", ".".join(map(str, sys.version_info[0:3]))
    print "Qt version:", QT_VERSION_STR
    print "VTK version:", vtkVersion.GetVTKVersion()

    import globvars as gv

    # configure program
    gv.configure()

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



if __name__ == '__main__':
    # -- initialize qt application here to avoid exit error
    app = QtGui.QApplication(sys.argv)
    main()
