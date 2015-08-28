' global variables of the program '
import platform
from datetime import datetime, timedelta
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
        # ======= updates:
        self._date_format = '%Y-%m-%d'
        self.updates_check = True
        self.updates_last = datetime.strptime('1900-1-1', self._date_format)
        self.updates_freq = timedelta(days=7)
        # ======= debugging:
        #save before each command execution
        self.debug_save_before = False
        self.debug_save_fn = '_debug.hmp'

        # ======= file paths:
        # Should be initialized from configure(fn)
        self.lib_crossgrid = "not-defined"
        self.opt_file = "not-defined"
        self.version = "not-defined"
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
        gridcom.BuildBoundaryGrid,
    ]

Flows = command.FlowCollection(_commands, framework.Framework)

#build window widgets: called after Flows building
#because it is used by widgets
mainWindow.initialize()


# =========================== Configuration
import xml.etree.ElementTree as ET


def configure():
    """ Reads configuration init from config_installed.py or config.py
        and fills prog_options and view_options.
        Checks for version updates if necessary.
    """
    import os.path
    import glob

    try:
        # If this is an installed version, then
        # current directory contains config_installed.py file
        import config_installed
        libdir = os.path.abspath(config_installed.libdir)
        homedir = os.path.abspath(config_installed.homedir)
        version = config_installed.version
    except:
        # For debug version
        import config
        libdir = os.path.abspath(config.libdir)
        homedir = os.path.abspath(config.homedir)
        version = config.version

    if not os.path.isdir(libdir):
        raise Exception('Home/Library directories are not correctly set')
    if not os.path.isdir(homedir):
        os.makedirs(homedir)

    # find crossgrid
    c = glob.glob("%s/*crossgrid*" % libdir)
    if len(c) == 0:
        raise Exception("libcrossgrid was not found at %s" % libdir)
    prog_options.lib_crossgrid = c[0]
    print '%s loaded' % c[0]

    #read options from Config.xml
    cnfn = "%s/Config.xml" % homedir
    try:
        xmlnode = ET.parse(cnfn).getroot()
        _read_prog_opt(xmlnode)
        _read_view_opt(xmlnode)
    except:
        print "HybMesh warning: Options were not read correctly"

    #rewrite options (if not all options were set from file)
    prog_options.opt_file = cnfn
    _rewrite_options()

    #updates
    prog_options.version = version
    _check_for_updates(version)

    #debug output directory
    prog_options.debug_save_fn = os.path.join(homedir,
            prog_options.debug_save_fn)
    if prog_options.debug_save_before:
        print ">>> Debug save to %s" % prog_options.debug_save_fn

    #change current directory to lib for Windows
    #for c libraries linkage
    if platform.system() == 'Windows':
        os.chdir(libdir)


def _read_prog_opt(xmlnode):
    'fills fields of prog_options struct'
    # debug
    n = xmlnode.find("PROGOPT/DEBUG")
    try:
        t = n.find("SAVE_BEFORE").text.strip()
        prog_options.debug_save_before = (t in ['True', 'true', '1'])
    except:
        print "HybMesh warning: failed reading DEBUG/SAVE_BEFORE"
    try:
        fn = n.find("SAVE_FN").text.strip()
        if len(fn) > 0:
            prog_options.debug_save_fn = fn
    except:
        print "HybMesh warning: failed reading DEBUG/SAVE_FN"

    # updates
    n = xmlnode.find('PROGOPT/UPDATES')
    try:
        t = n.find('CHECK').text.strip()
        prog_options.updates_check = (t in ['True', 'true', '1'])
    except:
        print "HybMesh warning: failed reading UPDATES/CHECK"
    try:
        t = int(n.find("FREQ").text.strip())
        prog_options.updates_freq = timedelta(days=t)
    except:
        print "HybMesh warning: failed reading UPDATES/FREQ"
    try:
        t = n.find("LAST").text.strip()
        t = datetime.strptime(t, prog_options._date_format)
        prog_options.updates_last = t
    except:
        print "HybMesh warning: failed reading UPDATES/LAST"


def _read_view_opt(xmlnode):
    'fills fields of view_options struct'
    n = xmlnode.find("VIEWOPT")
    try:
        t = map(int, n.find("GRID_COLOR").text.split())
        view_options.grid_color = (t[0], t[1], t[2])
    except:
        print "HybMesh warning: failed reading GRID_COLOR"
    try:
        t = map(int, n.find("CHECKED_COLOR").text.split())
        view_options.checked_color = (t[0], t[1], t[2])
    except:
        print "HybMesh warning: failed reading CHECKED_COLOR"
    try:
        t = int(n.find("GRID_LINE_WIDTH").text.strip())
        view_options.grid_line_width = t
    except:
        print "HybMesh warning: failed reading GRID_LINE_WIDTH"
    try:
        t = int(n.find("CONT_LINE_WIDTH").text.strip())
        view_options.cont_line_width = t
    except:
        print "HybMesh warning: failed reading CONT_LINE_WIDTH"
    try:
        t = int(n.find("CHECKED_GRID_LINE_WIDTH").text.strip())
        view_options.checked_grid_line_width = t
    except:
        print "HybMesh warning: failed reading CHECKED_CONT_LINE_WIDTH"
    try:
        t = int(n.find("CHECKED_CONT_LINE_WIDTH").text.strip())
        view_options.checked_cont_line_width = t
    except:
        pass


def _rewrite_options():
    'writes Config.xml from current prog/view_options'
    import bp
    outp = ET.Element("HybMeshConfig")
    _write_prog_opt(outp)
    _write_view_opt(outp)
    bp.xmlindent(outp)
    tree = ET.ElementTree(outp)
    tree.write(prog_options.opt_file, xml_declaration=True, encoding='utf-8')


def _write_prog_opt(xmlnode):
    'writes prog_options to xmlnode. Creates new nodes'
    n = ET.SubElement(xmlnode, "PROGOPT")
    d = prog_options
    #debug
    n1 = ET.SubElement(n, "DEBUG")
    ET.SubElement(n1, "SAVE_BEFORE").text = str(d.debug_save_before)
    ET.SubElement(n1, "SAVE_FN").text = str(d.debug_save_fn)
    #updates
    n1 = ET.SubElement(n, "UPDATES")
    ET.SubElement(n1, "CHECK").text = str(d.updates_check)
    ET.SubElement(n1, "LAST").text = d.updates_last.strftime(d._date_format)
    ET.SubElement(n1, "FREQ").text = str(d.updates_freq.days)


def _write_view_opt(xmlnode):
    'writes view_options to xmlnode. Creates new nodes'
    n = ET.SubElement(xmlnode, "VIEWOPT")
    d = view_options
    ET.SubElement(n, "GRID_COLOR").text = \
            ' '.join(map(str, d.grid_color))
    ET.SubElement(n, "CHECKED_COLOR").text = \
            ' '.join(map(str, d.checked_color))
    ET.SubElement(n, "GRID_LINE_WIDTH").text = str(d.grid_line_width)
    ET.SubElement(n, "CONT_LINE_WIDTH").text = str(d.cont_line_width)
    ET.SubElement(n, "CHECKED_GRID_LINE_WIDTH").text = \
            str(d.checked_grid_line_width)
    ET.SubElement(n, "CHECKED_CONT_LINE_WIDTH").text = \
            str(d.checked_cont_line_width)


def _get_last_version(rurl):
    import urllib2
    import json
    try:
        response = urllib2.urlopen(rurl, timeout=1)
        r = json.loads(response.read())
        return r['tag_name'], r['html_url']
    except:
        return (None, None)


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

    @classmethod
    def str_compare(cls, v1, v2):
        '(str, str) -> [-1 or 0 or 1]. Compares self and v2 string'
        [v1, v2] = [HybMeshVersion(v) for v in [v1, v2]]
        for i in range(min(len(v1.nums), len(v2.nums))):
            if v1.nums[i] < v2.nums[i]:
                return -1
            elif v1.nums[i] > v2.nums[i]:
                return 1
        return 0


def _compare_versions(cur_vers, (last_vers, url)):
    'if current version < last version raises message box'
    if last_vers is None:
        return
    cres = HybMeshVersion.str_compare(cur_vers, last_vers)
    if cres < 0:
        from PyQt4 import QtCore
        s = "New HybMesh release %s is available <a href='%s'>here</a>" \
            % (last_vers, url)
        mainWindow.emit(QtCore.SIGNAL(
            "inform(QString, QString)"), "Updates notification", s)


def _check_for_updates(cur_vers):
    #check updates options
    if not prog_options.updates_check:
        return
    #duration since last check
    d = datetime.now() - prog_options.updates_last
    if d < prog_options.updates_freq:
        return
    #make github request
    print "check for update"
    rurl = 'https://api.github.com/repos/' + \
        'kalininei/HybMesh/releases/latest'
    _compare_versions(cur_vers, _get_last_version(rurl))
    prog_options.updates_last = datetime.now()
    _rewrite_options()
