' global variables of the program '

import sys
from PyQt4.QtGui import QApplication
import mainwin
import command
import framework
import gridcom


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
        gridcom.RenameGrid2, gridcom.UniteGrids]

Flows = command.FlowCollection(_commands, framework.Framework)
