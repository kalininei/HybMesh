#!/usr/bin/env python
'import grid and contours'
import os.path
import copy
import ast
import itertools
import xml.etree.ElementTree as ET
from collections import OrderedDict
from PyQt4 import QtGui
import bp
import bgeom
import contour2
import grid2
import dlgs
import optview
import objcom


# ======================== Dialogs
class ImportGridDlg(dlgs.SimpleAbstractDialog):
    opts = OrderedDict([
            (0, 'native (*.hmg)'),
            (1, 'GridGen ASCII (*.net)')
        ])

    def __init__(self, parent=None):
        super(ImportGridDlg, self).__init__(parent)
        self.resize(400, 200)
        self.setWindowTitle("Import grid")

    def _default_odata(self, obj):
        "-> options struct with default values"
        obj.name = "ImportedGrid1"
        obj.filename = ""
        obj.tp = self.opts[0]

    def olist(self):
        "-> optview.OptionsList"
        return optview.OptionsList([
            ("Options", "Name",
                optview.SimpleOptionEntry(self.odata(), "name")),
            ("Format", "File type",
                optview.SingleChoiceOptionEntry(
                    self.odata(), "tp", self.opts.values())),
            ("Import from", "File",
                optview.OFileOptionEntry(self.odata(), "filename", None,
                    self._name_filters))
        ])

    def _name_filters(self):
        ret = []
        if self.odata().tp == self.opts[0]:
            ret.append("Native grid(*.hmg)")
        elif self.odata().tp == self.opts[1]:
            ret.append("GridGen grid(*.net)")
        ret.append("All files(*.*)")
        return ret

    def ret_value(self):
        '-> command.Command or None'
        od = copy.deepcopy(self.odata())
        #native
        if self.odata().tp == self.opts[0]:
            return ImportGridNative(od.filename, od.name)
        #*.net
        elif self.odata().tp == self.opts[1]:
            return ImportGridSimple34(od.filename, od.name)

    def check_input(self):
        "throws Exception if self.odata() has invalid fields"
        if not os.path.isfile(self.odata().filename):
            raise Exception("File doesn't exist")


class ImportContourDlg(dlgs.SimpleAbstractDialog):
    opts = OrderedDict([
            (0, 'native (*.hmc)'),
            (1, '"x y" ASCII'),
            (2, '"x y btype" ASCII'),
        ])

    def __init__(self, parent=None):
        super(ImportContourDlg, self).__init__(parent)
        self.resize(400, 300)
        self.setWindowTitle("Import contours")

    def _default_odata(self, obj):
        "-> options struct with default values"
        obj.name = "ImportedContour1"
        obj.filename = ""
        obj.tp = self.opts[1]
        obj.simplify = False
        obj.force_closed = True

    def olist(self):
        "-> optview.OptionsList"
        return optview.OptionsList([
            ("Options", "Name",
                optview.SimpleOptionEntry(self.odata(), "name")),
            ("Options", "Simplify",
                optview.BoolOptionEntry(self.odata(), "simplify")),
            ("Options", "Force closed",
                optview.BoolOptionEntry(self.odata(), "force_closed")),
            ("Format", "File type",
                optview.SingleChoiceOptionEntry(
                    self.odata(), "tp", self.opts.values())),
            ("Import from", "File",
                optview.OFileOptionEntry(self.odata(), "filename", None,
                    self._name_filters))
        ])

    def _name_filters(self):
        ret = []
        if self.odata().tp == self.opts[0]:
            ret.append("Native contours(*.hmc)")
        ret.append("All files(*.*)")
        return ret

    def ret_value(self):
        '-> command.Command or None'
        od = copy.deepcopy(self.odata())
        #native
        if self.odata().tp == self.opts[0]:
            return ImportContourNative(od.filename, od.name, od.simplify)
        #(x y)
        elif self.odata().tp == self.opts[1]:
            return ImportContourASCII(od.filename, od.name,
                    od.simplify, False, od.force_closed)
        #(x y b)
        elif self.odata().tp == self.opts[2]:
            return ImportContourASCII(od.filename, od.name,
                    od.simplify, True, od.force_closed)

    def check_input(self):
        "throws Exception if self.odata() has invalid fields"
        if not os.path.isfile(self.odata().filename):
            raise Exception("File doesn't exist")

    def _active_entries(self, entry):
        if entry.member_name in ["force_closed"]:
            return False
        else:
            return True


# ====================== Abstract Import Commands
class _AbstractImport(objcom.AbstractAddRemove):
    'Import basic class'
    def __init__(self, filename, name, argsdict):
        self.name = name
        self.filename = filename
        a = argsdict
        a['filename'] = filename
        a['name'] = name
        super(_AbstractImport, self).__init__(a)

    def _addrem_objects(self):
        try:
            self.__check_file_existance()
            cont, grid = self._read_data()
            if cont is not None:
                return [], [], [(self.name, cont)], []
            elif grid is not None:
                return [(self.name, grid)], [], [], []
        except Exception as e:
            import traceback
            print traceback.format_exc()
            QtGui.QMessageBox.critical(None, "Import failure",
                    str(e.message))
        return [], [], [], []

    def __check_file_existance(self):
        'checks self.filename existance and offers to enter new filename'
        if not os.path.isfile(self.filename):
            q = QtGui.QMessageBox.question(None, "File not found",
                    "No such file:\n%s.\n\nOpen another file?" % self.filename,
                    QtGui.QMessageBox.Open | QtGui.QMessageBox.Cancel)
            if q == QtGui.QMessageBox.Open:
                _, ext = os.path.splitext(self.filename)
                f = QtGui.QFileDialog.getOpenFileName(None,
                    'Load data', '', '(%s *%s);;All Files (*)' % (ext, ext))
                if f is not None:
                    self.filename = str(f)

        if not os.path.isfile(self.filename):
            raise Exception("File not found")

    #functions for overriding
    def _read_data(self):
        """ -> contour2.Contour2, grid2.Grid2.
            One of it should be None
        """
        raise NotImplementedError


class _AbstractImportContour(_AbstractImport):
    'Import contour basic class'
    def __init__(self, filename, name, simplify, argsdict):
        self.simplify = simplify
        a = argsdict
        a['simplify'] = simplify
        super(_AbstractImportContour, self).__init__(filename, name, a)

    def doc(self):
        return "Import contour from %s" % os.path.basename(self.filename)

    def _read_data(self):
        #build
        cont = self._read_contour()
        #simplify
        if self.simplify:
            c2 = cont.simplify(0)
            if c2 is not None:
                cont = c2
        return cont, None

    #functions for overriding
    def _read_contour(self):
        "-> contour2.Contour2 or raise an Exception"
        raise NotImplementedError


class _AbstractImportGrid(_AbstractImport):
    'Import contour basic class'
    def __init__(self, filename, name, argsdict={}):
        super(_AbstractImportGrid, self).__init__(filename, name, argsdict)

    def doc(self):
        return "Import grid from %s" % os.path.basename(self.filename)

    def _read_data(self):
        #build
        grid = self._read_grid()
        grid.build_contour()
        return None, grid

    @classmethod
    def fromstring(cls, slist):
        arg = ast.literal_eval(slist)
        filename, name = arg['filename'], arg['name']
        return cls(filename, name)

    #functions for overriding
    def _read_grid(self):
        "-> grid2.Grid2 or raise an Exception"
        raise NotImplementedError


#============================== Import contours
class ImportContourASCII(_AbstractImportContour):
    "plain ascii (x y) or (x y b)"
    def __init__(self, fn, nm, simp, readbnd=False, close=True):
        a = {'readbnd': readbnd, 'close': close}
        super(ImportContourASCII, self).__init__(fn, nm, simp, a)
        self.readbnd = readbnd
        self.close = close

    @classmethod
    def fromstring(cls, slist):
        arg = ast.literal_eval(slist)
        filename, name = arg['filename'], arg['name']
        simp = bp.dict_readbool(arg, 'simplify', False)
        rb = bp.dict_readbool(arg, 'readbnd', None)
        cl = bp.dict_readbool(arg, 'close', None)
        return cls(filename, name, simp, rb, cl)

    #overriding
    def _read_contour(self):
        st = map(float, file(self.filename).read().split())
        pts, eds, b = [], [], []
        it = iter(st)
        if self.readbnd:
            for x, y, t in zip(it, it, it):
                pts.append(bgeom.Point2(x, y))
                b.append(int(t))
        else:
            for x, y in zip(it, it):
                pts.append(bgeom.Point2(x, y))
                b.append(0)

        if len(pts) < 2:
            raise Exception("Too few points defined")

        closed = self.close
        if pts[-1].x == pts[0].x and pts[-1].y == pts[0].y:
            closed = True
            pts = pts[:-1]
            b = b[:-1]

        if closed and len(pts) < 3:
            raise Exception("Too few points defined")

        eds = [[i, i + 1] for i in range(len(pts))]
        if closed:
            eds[-1][1] = 0
        else:
            eds = eds[:-1]

        #return
        return contour2.Contour2.create_from_point_set(pts, eds, b)


class ImportContourNative(_AbstractImportContour):
    "import from *.hmc"
    def __init__(self, fn, nm, simp):
        super(ImportContourNative, self).__init__(fn, nm, simp, {})

    @classmethod
    def fromstring(cls, slist):
        arg = ast.literal_eval(slist)
        filename, name = arg['filename'], arg['name']
        simp = bp.dict_readbool(arg, 'simplify', False)
        return cls(filename, name, simp)

    #overriding
    def _read_contour(self):
        xmlnode = ET.parse(self.filename).getroot()
        return contour2.Contour2.create_from_xml(xmlnode)


# ======================= Import grids
class ImportGridNative(_AbstractImportGrid):
    "from native hmg format"
    def __init__(self, fn, nm):
        super(ImportGridNative, self).__init__(fn, nm)

    #overriding
    def _read_grid(self):
        xmlnode = ET.parse(self.filename).getroot()
        return grid2.Grid2.create_from_xml(xmlnode)


class ImportGridSimple34(_AbstractImportGrid):
    "from *.net format"
    def __init__(self, fn, nm):
        super(ImportGridSimple34, self).__init__(fn, nm)

    #overriding
    def _read_grid(self):
        dt = file(self.filename).readlines()
        pts, cells3, cells4, pbnd = [], [], [], {}
        it = iter(dt)

        #numbers
        for v in it:
            d = v.split()
            if len(d) >= 3:
                [npt, nc3, nc4] = map(int, d[:3])
                break

        #points
        while len(pts) < npt:
            d = next(it).split()
            if len(d) >= 4:
                pts.append(bgeom.Point2(float(d[0]), float(d[1])))
                b = int(d[2]) - 1   # using b-1 since all btypes start with 1
                if b > 0:
                    pbnd[len(pts) - 1] = b

        #cells3
        while len(cells3) < nc3:
            d = next(it).split()
            if len(d) >= 3:
                cells3.append(map(lambda s: int(s) - 1, d[:3]))

        #cells4
        while len(cells4) < nc4:
            d = next(it).split()
            if len(d) >= 4:
                #if Six or Ten - elements are given -> ignore
                if len(d) > 5:
                    break
                cells4.append(map(lambda s: int(s) - 1, d[:4]))

        cells = cells4 if len(cells4) > 0 else cells3
        g = grid2.Grid2.from_points_cells2(pts, cells)

        #boundaries
        bc = g.boundary_contours()
        bset = {}
        for eind in itertools.chain.from_iterable(bc):
            e = g.edges[eind]
            if e[0] in pbnd:
                bset[eind] = pbnd[e[0]]
            elif e[1] in pbnd:
                bset[eind] = pbnd[e[1]]
        g.set_edge_bnd(bset)
        return g
