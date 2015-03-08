' global variables of the program '

import sys
from PyQt4.QtGui import QApplication
import command
import framework
import gridcom
import mainwin


def actual_flow():
    ' -> command.CommandFlow. Returns actual flow '
    return Flows.get_actual_flow()


def actual_data():
    ' ->framework.Framework. Returns actual data '
    return Flows.get_actual_flow()._receiver

# -- initialize qt application
app = QApplication(sys.argv)

#build main window
mainWindow = mainwin.MainWindow()

#build flow collection
#it is build after mainwindow because FrameworkVis should use mainWindow
_commands = [gridcom.AddUnfRectGrid, gridcom.RemoveGrid2,
        gridcom.AddUnfCircGrid, gridcom.MoveGrids, gridcom.RotateGrids,
        gridcom.RenameGrid2, gridcom.UniteGrids, gridcom.AddUnfRingGrid,
        gridcom.ScaleGrids]

Flows = command.FlowCollection(_commands, framework.Framework)

#build window widgets: called after Flows building
#because it is used by widgets
mainWindow.initialize()
