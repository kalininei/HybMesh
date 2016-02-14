import os.path
import itertools
from HybMeshPyPack import basic, gdata
import HybMeshPyPack.basic.geom as bg
from HybMeshPyPack.gdata import contour2
from HybMeshPyPack.gdata import grid2
from HybMeshPyPack.basic.geom import Point2
import contcom
import command
import objcom
import re


# ====================== Abstract Import Commands
class _AbstractImport(objcom.AbstractAddRemove):
    'Import basic class'
    def __init__(self, argsdict):
        if 'name' not in argsdict:
            argsdict['name'] = argsdict['filename']
        super(_AbstractImport, self).__init__(argsdict)

    def _addrem_objects(self):
        self.__check_file_existance()
        try:
            cont, grid = self._read_data()
        except Exception as e:
            raise command.ExecutionError("Import failure", self, e)
        if cont is not None:
            return [], [], [(self.options['name'], cont)], []
        elif grid is not None:
            return [(self.options['name'], grid)], [], [], []

    def __check_file_existance(self):
        'checks self.filename existance and offers to enter new filename'
        fn = self.options['filename']
        if not os.path.isfile(fn):
            raise command.ExecutionError("File not found %s" % fn, self)

    #functions for overriding
    def _read_data(self):
        """ -> contour2.Contour2, grid2.Grid2.
            One of it should be None
        """
        raise NotImplementedError


class _AbstractImportContour(_AbstractImport):
    'Import contour basic class'
    def __init__(self, argsdict):
        if 'simplify' not in argsdict:
            argsdict['simplify'] = False
        super(_AbstractImportContour, self).__init__(argsdict)

    def doc(self):
        return "Import contour from %s" % os.path.basename(
                self.options['filename'])

    def _read_data(self):
        #build
        cont = self._read_contour()
        #simplify
        if self.options['simplify']:
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
    def __init__(self, argsdict):
        super(_AbstractImportGrid, self).__init__(argsdict)

    def doc(self):
        return "Import grid from %s" % os.path.basename(
                self.options['filename'])

    def _read_data(self):
        #build
        grid = self._read_grid()
        grid.build_contour()
        return None, grid

    #functions for overriding
    def _read_grid(self):
        "-> grid2.Grid2 or raise an Exception"
        raise NotImplementedError


# ========================== Import contours
class ImportContourASCII(_AbstractImportContour):
    "plain ascii (x y) or (x y b)"
    def __init__(self, arg):
        if 'btypes' not in arg:
            arg['btypes'] = False
        if 'force_closed' not in arg:
            arg['force_closed'] = False
        super(ImportContourASCII, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        """ name - new contour name,
            filename - filename
            simplify - bool
            btypes - if false read (x y). (x y boundary_type_integer) otherwise
            force_closed - treat contour as closed even if end points are not
                    equal
        """
        return {'name': command.BasicOption(str),
                'filename': command.BasicOption(str),
                'simplify': command.BoolOption(),
                'btypes': command.BoolOption(),
                'force_closed': command.BoolOption()
                }

    def _read_contour(self):
        nm, fn, bt, fc = [self.options[x] for x in
                ['name', 'filename', 'btypes', 'force_closed']]
        st = map(float, file(fn).read().split())
        pts, eds, b = [], [], []
        it = iter(st)
        if bt:
            for x, y, t in zip(it, it, it):
                pts.append(bg.Point2(x, y))
                b.append(int(t))
        else:
            for x, y in zip(it, it):
                pts.append(bg.Point2(x, y))
                b.append(0)

        if len(pts) < 2:
            raise command.ExecutionError("Too few points defined", self)

        closed = fc
        if pts[-1].x == pts[0].x and pts[-1].y == pts[0].y:
            closed = True
            pts = pts[:-1]
            b = b[:-1]

        if closed and len(pts) < 3:
            raise command.ExecutionError("Too few points defined", self)

        eds = [[i, i + 1] for i in range(len(pts))]
        if closed:
            eds[-1][1] = 0
        else:
            eds = eds[:-1]

        #return
        return gdata.contour2.Contour2.create_from_point_set(pts, eds, b)


class ImportContourNative(_AbstractImportContour):
    "import from *.hmc"
    def __init__(self, arg):
        super(ImportContourNative, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        """ name - new contour name,
            filename - filename
            simplify - bool
        """
        return {'name': command.BasicOption(str),
                'filename': command.BasicOption(str),
                'simplify': command.BoolOption(),
                }

    #overriding
    def _read_contour(self):
        import xml.etree.ElementTree as ET
        xmlnode = ET.parse(self.options['filename']).getroot()
        return gdata.contour2.Contour2.create_from_xml(xmlnode)


#========================= Import grids
class ImportGridNative(_AbstractImportGrid):
    "from native hmg format"
    def __init__(self, arg):
        super(ImportGridNative, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        """ name - new grid name,
            filename - filename
        """
        return {'name': command.BasicOption(str),
                'filename': command.BasicOption(str),
                }

    #overriding
    def _read_grid(self):
        import xml.etree.ElementTree as ET
        xmlnode = ET.parse(self.options['filename']).getroot()
        return gdata.grid2.Grid2.create_from_xml(xmlnode)


class ImportGridMSH(_AbstractImportGrid):
    "from *.msh format"
    def __init__(self, arg):
        super(ImportGridMSH, self).__init__(arg)
        # flow commands which were executed during grid reading
        self.addcom = []

    @classmethod
    def _arguments_types(cls):
        """ name - new grid name,
            filename - filename
        """
        return {'name': command.BasicOption(str),
                'filename': command.BasicOption(str),
                }

    @staticmethod
    def _parse_parant(txt):
        """ returns list of strings within first level paranthesis """
        ret = []
        parant = {}
        f = -1
        while 1:
            f = txt.find('(', f + 1)
            if f == -1:
                break
            parant[f] = '('
        while 1:
            f = txt.find(')', f + 1)
            if f == -1:
                break
            parant[f] = ')'

        oparant = []
        for k in sorted(parant.keys()):
            oparant.append((k, parant[k]))
        st = -1
        counter = 0
        for p in oparant:
            if st == -1 and p[1] == '(':
                st = p[0]
                counter = 0
            else:
                if p[1] == '(':
                    counter = counter + 1
                elif p[1] == ')':
                    if counter > 0:
                        counter = counter - 1
                    else:
                        ret.append(txt[st + 1:p[0]])
                        st = -1
        return ret

    @staticmethod
    def _parse_hex(txt):
        """ -> [list of integers] from text """
        return [int(x, 16) for x in txt.split()]

    @staticmethod
    def _parse_floats(txt):
        """ -> [list of floats] """
        return map(float, txt.split())

    def _find_highest_bndindex(self):
        bset = self.receiver.get_bnd_types()._ind_set()
        return max(bset)

    #overriding
    def _read_grid(self):
        fn = self.options['filename']
        dt = file(fn).read()
        dt = self._parse_parant(dt)
        dtp = {}
        for d in dt:
            k = re.search('\d+', d).end(0)
            w1, w2 = d[:k], d[k+1:]
            w1 = int(w1)
            if w1 != 0:
                if w1 not in dtp:
                    dtp[w1] = []
                dtp[w1].append(w2)
        # check dimension
        if (int(dtp[2][0]) != 2):
            raise "Invalid msh file dimension"

        # read nodes
        points = []
        for plines in dtp[10]:
            ln = self._parse_parant(plines)
            if len(ln) < 2:
                continue
            info = self._parse_hex(ln[0])
            index0, index1 = info[1], info[2]
            if (len(points) < index1):
                points.extend([None] * (index1 - len(points)))
            coords = self._parse_floats(ln[1])
            for i in range(index1 - index0 + 1):
                p = Point2(coords[2 * i], coords[2 * i + 1])
                points[i + index0 - 1] = p

        # read edges as [p0, p1, cell_right, cell_left, zone_type]
        eds = []
        for clines in dtp[13]:
            ln = self._parse_parant(clines)
            info = self._parse_hex(ln[0])
            if info[0] == 0:
                continue
            index0, index1 = info[1], info[2]
            if (len(eds) < index1):
                eds.extend([None] * (index1 - len(eds)))
            conn = self._parse_hex(ln[1])
            for i in range(index1 - index0 + 1):
                if info[4] == 2:
                    line = conn[4 * i:4 * i + 4]
                elif info[4] == 0:
                    line = conn[5 * i + 1: 5 * i + 5]
                eds[i + index0 - 1] = [line[0] - 1,
                        line[1] - 1,
                        line[2] - 1,
                        line[3] - 1,
                        info[0]]

        # read boundary names
        indshift = self._find_highest_bndindex() + 1
        ztps = {}
        for zline in dtp[45]:
            w = self._parse_parant(zline)[0].split()
            if w[1] != "interior":
                ztps[int(w[0])] = (w[2], indshift)
                indshift += 1
        # adjust b types to unique ones with respect to workflow
        # filter out bc which are not present
        ztps2 = {}
        for e in eds:
            if e[4] in ztps:
                ztps2[ztps[e[4]][1]] = ztps[e[4]][0]
                e[4] = ztps[e[4]][1]
        ztps = ztps2

        # insert boundary types to the workflow
        for k, v in ztps.iteritems():
            self.addcom.append(
                    contcom.EditBoundaryType({"index": k, "name": v}))
            self.addcom[-1].do(self.receiver)

        # create grid
        return gdata.grid2.Grid2.from_points_edges(points, eds)

    # overriding since additional commands present in reading grid procedure
    def _clear(self):
        super(ImportGridMSH, self)._clear()
        self.addcom = []

    def _undo(self):
        super(ImportGridMSH, self)._undo()
        for c in self.addcom:
            c._undo()

    def _redo(self):
        super(ImportGridMSH, self)._redo()
        for c in self.addcom:
            c._redo()


class ImportGridGMSH(_AbstractImportGrid):
    "from gmsh format"
    def __init__(self, arg):
        super(ImportGridGMSH, self).__init__(arg)
        # flow commands which were executed during grid reading
        self.addcom = []

    @classmethod
    def _arguments_types(cls):
        """ name - new grid name,
            filename - filename
        """
        return {'name': command.BasicOption(str),
                'filename': command.BasicOption(str),
                }

    def _find_highest_bndindex(self):
        bset = self.receiver.get_bnd_types()._ind_set()
        return max(bset)

    def _parse(self):
        """ ->
            {bindex: bname}
            [ ... points ... ],
            [ ... 2d cells ... ],
            [ [p1, p2, btype], ... ]
        """
        r1, r2, r3, r4 = {}, [], [], []
        r2m = {}
        fn = self.options['filename']
        dt = file(fn).readlines()
        dt = map(lambda x: x.strip(), dt)
        dt = [item for item in dt if len(item) > 0]
        i = 0
        while i < len(dt):
            # 1
            if dt[i] == "$PhysicalNames":
                k = int(dt[i + 1])
                for d in dt[i + 2:i + 2 + k]:
                    ds = d.split(' ', 2)
                    r1[int(ds[1])] = ds[2].strip('"')
                i += k + 1
            # 2
            if dt[i] == "$Nodes":
                k = int(dt[i + 1])
                for d in dt[i + 2:i + 2 + k]:
                    ds = d.split()
                    r2m[int(ds[0])] = Point2(float(ds[1]), float(ds[2]))
                i += k + 1
            # 3
            if dt[i] == "$Elements":
                k = int(dt[i + 1])
                for d in dt[i + 2:i + 2 + k]:
                    ds = map(int, d.split())
                    if ds[1] == 1:
                        # edge
                        tp = ds[3] if ds[2] > 0 else 0
                        r4.append([ds[-2], ds[-1], tp])
                    elif ds[1] == 2:
                        # triangle cell
                        r3.append([ds[-3], ds[-2], ds[-1]])
                    elif ds[1] == 3:
                        # quad cell
                        r3.append([ds[-4], ds[-3], ds[-2], ds[-1]])
                    elif ds[1] == 15:
                        # ignore point cell
                        pass
                    else:
                        raise Exception("Not supported gmsh element type")
                i += k + 1
            i += 1
        # point renumbering for successive indexing starting from zero
        oldnew = {}
        i = 0
        for ind in sorted(r2m.keys()):
            oldnew[ind] = i
            r2.append(r2m[ind])
            i += 1
        for a in r4:
            a[0] = oldnew[a[0]]
            a[1] = oldnew[a[1]]
        for i, a in enumerate(r3):
            r3[i] = map(lambda x: oldnew[x], a)
        return r1, r2, r3, r4

    def _build_bmap(self, bt, bedges):
        """ -> { gmsh_type: (new_type, new_name) } """
        bset = set()
        for b in bedges:
            bset.add(b[2])
        bshift = self._find_highest_bndindex()
        ret = {}
        i = 1
        for k in sorted(bt.keys()):
            if k in bset:
                ret[k] = (bshift + i, bt[k])
                i += 1
        for k in bset:
            if k not in ret:
                ret[k] = (bshift + i, "Unknown")
                i += 1
        return ret

    def _add_boundaries(self, bmap):
        for k, v in bmap.iteritems():
            self.addcom.append(
                    contcom.EditBoundaryType({
                        "index": v[0], "name": v[1]}))
            self.addcom[-1].do(self.receiver)

    @staticmethod
    def _set_boundaries(grid, bmap, bedges):
        # edge index -> boundary index
        imap = {}
        bed1 = grid.boundary_contours()
        bed = []
        for b in bed1:
            bed.extend(b)
        for i, b in enumerate(bed):
            for b2 in bedges:
                if b2[0] in grid.edges[b] and b2[1] in grid.edges[b]:
                    imap[i] = bmap[b2[2]][0]
        # set
        grid.set_edge_bnd(imap)

    #overriding
    def _read_grid(self):
        bt, pts, cls, bedges = self._parse()
        ret = gdata.grid2.Grid2.from_points_cells2(pts, cls)
        bmap = self._build_bmap(bt, bedges)
        self._add_boundaries(bmap)
        self._set_boundaries(ret, bmap, bedges)
        ret.build_contour()
        return ret

    # overriding since additional commands present in reading grid procedure
    def _clear(self):
        super(ImportGridGMSH, self)._clear()
        self.addcom = []

    def _undo(self):
        super(ImportGridGMSH, self)._undo()
        for c in self.addcom:
            c._undo()

    def _redo(self):
        super(ImportGridGMSH, self)._redo()
        for c in self.addcom:
            c._redo()


class ImportGridSimple34(_AbstractImportGrid):
    "from *.net gridgen format"
    def __init__(self, arg):
        super(ImportGridSimple34, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        """ name - new grid name,
            filename - filename
        """
        return {'name': command.BasicOption(str),
                'filename': command.BasicOption(str),
                }

    #overriding
    def _read_grid(self):
        fn = self.options['filename']
        dt = file(fn).readlines()
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
                pts.append(bg.Point2(float(d[0]), float(d[1])))
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
        g = gdata.grid2.Grid2.from_points_cells2(pts, cells)

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
