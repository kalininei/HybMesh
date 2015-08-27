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
            (3, 'Fluent mesh (*.msh)'),
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
        elif self.odata().tp == self.opts[3]:
            ret.append("msh files(*.msh)")
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
        #msh format
        elif self.odata().tp == self.opts[3]:
            ex = ExportGridFluent(od.filename, g)

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


class ExportGridFluent(AbstractExportGrid):
    'Export grid to *.msh format'
    def __init__(self, filename, grid):
        super(ExportGridFluent, self).__init__(filename, grid)

    def export(self):
        import globvars

        def toh(i):
            '->str. integer to hex format'
            #return str(i)
            return hex(i)[2:]

        #check if all cells are tri/tetra
        cellsinfo = self.grid.cell_types_info()
        c3 = cellsinfo[3] if 3 in cellsinfo else 0
        c4 = cellsinfo[4] if 4 in cellsinfo else 0
        aa = self.grid.cells_nodes_connect()
        for a in aa:
            if len(a) == 5:
                print "Start"
                for x in a:
                    print " %i: %f %f" % (x, self.grid.points[x].x, self.grid.points[x].y)
        if c3 + c4 != self.grid.n_cells():
            raise Exception("Only triangle/tetrahedral grids could be" +
                    " exported to Fluent *.msh format")
        out = [
                """(0 "HybMesh to Fluent File")""",
                """(0 "Dimensions")""",
                """(2 2)""",
        ]
        #Zones:
        #0    - no zone
        #1    - vericies default
        #2    - fluid for cells
        #3    - default interior
        #4..N - bc's

        #faces by zone
        ns = {i: 0 if (c[0] == -1 or c[1] == -1) else -1
                for i, c in enumerate(self.grid.edges_cells_connect())}
        for k, v in self.grid.bt.iteritems():
            ns[k] = v
        tps, tpind = [], [-1]
        for v in ns.values():
            if v not in tpind:
                tpind.append(v)
        for v in tpind:
            tps.append([x[0] for x in filter(lambda x: x[1] == v,
                ns.iteritems())])

        #verticies
        out.extend([
                """(0 "Verticies")""",
                "(10 (0 1 %s 1 2))" % toh(self.grid.n_points()),
                "(10 (1 1 %s 1 2)(" % toh(self.grid.n_points()),
        ])
        out.extend(['  %16.12e %16.12e' % (p.x, p.y)
            for p in self.grid.points])
        out.extend(["))"])
        #faces
        fconnect = self.grid.edges_cells_connect()
        out.extend([
            """(0 "Faces")""",
            "(13 (0 1 %s 0))" % toh(self.grid.n_edges()),
        ])
        #interior faces
        out.append("(13 (3 1 %s 2 0)(" % toh(len(tps[0])))
        out.extend(["2 %s %s %s %s" % (
                toh(self.grid.edges[i][0] + 1),
                toh(self.grid.edges[i][1] + 1),
                toh(fconnect[i][0] + 1),
                toh(fconnect[i][1] + 1)
            ) for i in tps[0]])
        out.append('))')
        #boundary nodes
        c0 = len(tps[0])
        for i, t in enumerate(tps[1:]):
            c1 = c0 + len(t)
            out.append("(13 (%s %s %s 3 0)(" % (
                         toh(4 + i), toh(c0+1), toh(c1)
                    )
            )
            out.extend(["2 %s %s %s %s" % (
                        toh(self.grid.edges[i][0] + 1),
                        toh(self.grid.edges[i][1] + 1),
                        toh(fconnect[i][0] + 1),
                        toh(fconnect[i][1] + 1))
                for i in t])
            out.append('))')
            c0 = c1

        #cells
        out.extend([
            """(0 "Cells")""",
            """(12 (0 1 %s 0))""" % toh(self.grid.n_cells()),
            """(12 (2 1 %s 1 0)(""" % toh(self.grid.n_cells()),
        ])
        out.extend([('1' if len(c) == 3 else '3') for c in self.grid.cells])
        out.append('))')
        #zones
        out.extend([
            """(0 "Zones")""",
            "(45 (2 fluid fluid)())",
            "(45 (3 interior default-interior)())"
        ])
        out.extend([
            "(45 (%s wall %s)())" % (
                toh(4 + i),
                globvars.actual_data().boundary_types.get(index=ti).name
            )
            for i, ti in enumerate(tpind[1:])])
        self._write_to_textfile(out)
