' global variables of the program '
import platform
import command
import framework
import gridcom
import objcom
import contcom
import importcom
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
        self.debug_save_before = True
        self.debug_save_fn = '_debug.hmp'

        # ======= external libraries:
        # Should be initialized from configure
        self.lib_crossgrid = "not-defined"
prog_options = ProgOptions()

#build main window
mainWindow = mainwin.MainWindow()

#build flow collection
#it is build after mainwindow because FrameworkVis should use mainWindow
_commands = [
        gridcom.AddUnfRectGrid,
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
        contcom.SetBTypeToContour,
        importcom.ImportContourNative,
        importcom.ImportContourASCII,
        importcom.ImportGridNative,
        importcom.ImportGridSimple34,
    ]

Flows = command.FlowCollection(_commands, framework.Framework)

#build window widgets: called after Flows building
#because it is used by widgets
mainWindow.initialize()


def configure(fn):
    """ reads configure file and
        fils ProgOptions and ViewOptions
    """
    import os.path
    import glob
    homedir, libdir = None, None
    confs = file(fn, 'r').readlines()
    for s in confs:
        k = s.strip().split(None, 1)
        if len(k) < 2 or k[0][0] == '#':
            continue
        elif k[0] == 'libdir':
            libdir = os.path.expanduser(k[1])
        elif k[0] == 'homedir':
            homedir = os.path.expanduser(k[1])
    if homedir is None or libdir is None or not os.path.isdir(libdir):
        raise Exception('Home/Library directories are not correctly set')
    if not os.path.isdir(homedir):
        os.makedirs(homedir)

    # find crossgrid
    c = glob.glob("%s/*crossgrid*" % libdir)
    if len(c) == 0:
        raise Exception("libcrossgrid was not found at %s" % libdir)
    prog_options.lib_crossgrid = c[0]
    print '%s loaded' % c[0]

    #debug output directory
    prog_options.debug_save_fn = os.path.join(homedir,
            prog_options.debug_save_fn)
    
    #change current directory to lib for Windows
    #for c libraries linkage
    if platform.system() == 'Windows':
        os.chdir(libdir)
