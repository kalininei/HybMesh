' global variables of the program '

import command
import framework
import gridcom
import objcom
import contcom
import mainwin


def actual_flow():
    ' -> command.CommandFlow. Returns actual flow '
    return Flows.get_actual_flow()


def actual_data():
    ' ->framework.Framework. Returns actual data '
    return Flows.get_actual_flow()._receiver


def mainvtk_render():
    "calls Render() command to mainWindoe.vtkWidget"
    mainWindow.vtkWidget.GetRenderWindow().Render()


# ==================== Program options
class ViewOptions(object):
    'View Options struct'
    def __init__(self):
        #colors
        self.checked_color = (255, 0, 0)
        self.grid_color = (180, 180, 180)
        #line widths
        self.grid_line_width = 1
        self.cont_line_width = 1
        self.checked_grid_line_width = 2
        self.checked_cont_line_width = 3
view_options = ViewOptions()


class ProgOptions(object):
    'Options of the program behavior'
    def __init__(self):
        # ======= debugging:
        #save before each command execution
        self.debug_save_before = False
        self.debug_save_fn = '../_debug.cgp'
prog_options = ProgOptions()

#build main window
mainWindow = mainwin.MainWindow()

#build flow collection
#it is build after mainwindow because FrameworkVis should use mainWindow
_commands = [gridcom.AddUnfRectGrid,
        gridcom.AddUnfCircGrid,
        gridcom.AddUnfRingGrid,
        gridcom.UniteGrids,
        gridcom.RenameGrid,
        gridcom.ExcludeContours,
        objcom.RemoveGrid,
        objcom.MoveGrids,
        objcom.RotateGrids,
        objcom.ScaleGrids,
        objcom.CopyGrid,
        contcom.RenameContour,
        contcom.AddRectCont,
        contcom.UniteContours,
        contcom.EditBoundaryType,
        contcom.GridBndToContour,
        contcom.SimplifyContours,
        contcom.SetBTypeToContour]

Flows = command.FlowCollection(_commands, framework.Framework)

#build window widgets: called after Flows building
#because it is used by widgets
mainWindow.initialize()
