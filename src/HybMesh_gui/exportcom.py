#!/usr/bin/env python
'export grid and contours'
import copy
import itertools
from collections import OrderedDict
import xml.etree.ElementTree as ET
from PyQt4 import QtGui
import bp
import dlgs
import optview


class ExportGridDlg(dlgs.SimpleAbstractDialog):
    opts = OrderedDict([
            (0, 'native (*.hmg)'),
            (2, 'vtk ASCII (*.vtk)'),
            (1, 'GridGen ASCII (*.net)'),
        ])

    def __init__(self, grd, all_grids, parent=None):
        """ grd - name of exporting grid or None
            all_grids - name->grid2 set
        """
        self.all_grids = all_grids
        self.odata().grd = grd if grd is not None else\
                all_grids.keys()[0] if len(all_grids) > 0 else ""
        super(ExportGridDlg, self).__init__(parent)
        self.resize(400, 200)
        self.setWindowTitle("Export grid")

    def _default_odata(self, obj):
        "-> options struct with default values"
        obj.name = "ImportedGrid1"
        obj.filename = ""
        obj.tp = self.opts[0]
        obj.grd = ""

    def olist(self):
        "-> optview.OptionsList"
        return optview.OptionsList([
            ("Source", "Grid",
                optview.SingleChoiceOptionEntry(
                    self.odata(), "grd", self.all_grids.keys())),
            ("Format", "File type",
                optview.SingleChoiceOptionEntry(
                    self.odata(), "tp", self.opts.values())),
            ("Export to", "File",
                optview.SFileOptionEntry(self.odata(), "filename",
                    filter_func=self._name_filters,
                    defname_func=self._get_gname)),
        ])

    def _name_filters(self):
        ret = []
        if self.odata().tp == self.opts[0]:
            ret.append("Native grid(*.hmg)")
        elif self.odata().tp == self.opts[1]:
            ret.append("GridGen grid(*.net)")
        elif self.odata().tp == self.opts[2]:
            ret.append("vtk files(*.vtk)")
        ret.append("All files(*.*)")
        return ret

    def _get_gname(self):
        return self.odata().grd

    def ret_value(self):
        '-> None. Saves to file'
        od = copy.deepcopy(self.odata())
        g = self.all_grids[od.grd]
        ex = None
        #native
        if self.odata().tp == self.opts[0]:
            ex = ExportGridNative(od.filename, g)
        #GridGen
        elif self.odata().tp == self.opts[1]:
            ex = ExportGridGridGen(od.filename, g)
        #vtk ascii
        elif self.odata().tp == self.opts[2]:
            ex = ExportGridVTK(od.filename, g)

        if ex is not None:
            try:
                ex.export()
            except Exception as e:
                import traceback
                print traceback.format_exc()
                QtGui.QMessageBox.critical(None, "Export failure",
                        str(e))



# ====================== Abstract export
class AbstractExportGrid(object):
    'Export grid to hmg format'
    def __init__(self, filename, grid):
        self.filename = filename
        self.grid = grid

    #function for override
    def export(self):
        """ writes grid to file or generates exception
        """
        raise NotImplementedError

    def _write_to_textfile(self, out):
        'writes list of string to file'
        #write to file
        out = [s + '\n' for s in out]
        outf = file(self.filename, 'w')
        outf.writelines(out)
        outf.close()


# ====================== Export grid
class ExportGridNative(AbstractExportGrid):
    'Export grid to hmg format'
    def __init__(self, filename, grid):
        super(ExportGridNative, self).__init__(filename, grid)

    def export(self):
        """ writes grid to file or generates exception
        """
        outp = ET.Element("GRID2")
        self.grid.xml_save(outp)
        bp.xmlindent(outp)
        tree = ET.ElementTree(outp)
        tree.write(self.filename, xml_declaration=True, encoding='utf-8')


class ExportGridGridGen(AbstractExportGrid):
    'Export grid to *.net format'
    def __init__(self, filename, grid):
        super(ExportGridGridGen, self).__init__(filename, grid)

    def export(self):
        """ writes grid to file or generates exception
        """
        cellsinfo = self.grid.cell_types_info()
        c3 = cellsinfo[3] if 3 in cellsinfo else 0
        c4 = cellsinfo[4] if 4 in cellsinfo else 0
        if c3 + c4 != self.grid.n_cells():
            raise Exception("Only triangle/tetrahedral grids could be" +
                    " exported to *.net format")
        #points boundaries
        bnd = [0] * self.grid.n_points()
        for e in itertools.chain.from_iterable(self.grid.boundary_contours()):
            bnd[self.grid.edges[e][0]] = 1
            bnd[self.grid.edges[e][1]] = 1

        for e, b in self.grid.bt.items():
            bnd[self.grid.edges[e][0]] = b + 1
            bnd[self.grid.edges[e][1]] = b + 1

        out = []
        cn = self.grid.cells_nodes_connect()

        def write_points(pts):
            out.extend(
                    ['%s %i %i' % (str(p), b, i + 1) for i, (p, b) in
                        enumerate(zip(pts, bnd))]
                )

        def write_eds3(eds):
            out.extend(
                    ['%i %i %i  %i' % (c[0] + 1, c[1] + 1, c[2] + 1, i + 1)
                        for i, c in enumerate(eds)]
                )

        if c4 == 0:
            out.append('%i %i %i' % (self.grid.n_points(), c3, 0))
            write_points(self.grid.points)
            write_eds3(cn)
        else:
            cells3 = []
            if c3 == 0:
                for c in cn:
                    cells3.extend([[c[0], c[1], c[2]], [c[0], c[2], c[3]]])
                out.append('%i %i %i' % (self.grid.n_points(), 2 * c4, c4))
                write_points(self.grid.points)
                write_eds3(cells3)
                out.extend(
                        ['%i %i %i %i  %i' % (c[0] + 1, c[1] + 1,
                            c[2] + 1, c[3] + 1, i + 1)
                            for i, c in enumerate(cn)]
                    )
            else:
                for c in cn:
                    if len(c) == 3:
                        cells3.append(c)
                    else:
                        cells3.extend([[c[0], c[1], c[2]], [c[0], c[2], c[3]]])
                out.append('%i %i %i' % (self.grid.n_points(), len(cells3), 0))
                write_points(self.grid.points)
                write_eds3(cells3)

        self._write_to_textfile(out)


class ExportGridVTK(AbstractExportGrid):
    'Export grid to *.vtk format'
    def __init__(self, filename, grid):
        super(ExportGridVTK, self).__init__(filename, grid)

    def export(self):
        out = [
            "# vtk DataFile Version 3.0",
            self.grid.__class__.__name__, "ASCII"]
        #Points
        out.append("DATASET UNSTRUCTURED_GRID")
        out.append("POINTS %i float" % self.grid.n_points())
        out.extend(
                ['%s 0' % str(p) for p in self.grid.points]
            )
        #Cells
        cp = self.grid.cells_nodes_connect()
        cd = sum(len(c) for c in cp) + len(cp)
        out.append("CELLS %i %i" % (len(cp), cd))
        out.extend(
                ['%i %s' % (len(c), ' '.join(map(str, c)))
                    for c in cp]
            )
        #CellTypes
        out.append("CELL_TYPES %i" % self.grid.n_cells())
        out.extend(['7'] * self.grid.n_cells())

        self._write_to_textfile(out)

