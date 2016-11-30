import os.path
import re
from hybmeshpack import gdata
import hybmeshpack.basic.geom as bg
from hybmeshpack.basic.geom import Point2
import contcom
import command
import objcom
from hybmeshpack import hmcore as hmcore
from hybmeshpack.hmcore import g2 as g2core
from hybmeshpack.hmcore import c2 as c2core
from hybmeshpack.hmcore import g3 as g3core
from hybmeshpack.hmcore import s3 as s3core


# ====================== Abstract Import Commands
class _AbstractImport(objcom.AbstractAddRemove):
    'Import basic class'
    def __init__(self, argsdict):
        if 'name' not in argsdict:
            argsdict['name'] = ""
        super(_AbstractImport, self).__init__(argsdict)

    def _addrem_objects(self):
        self.__check_file_existance()
        g2list, c2list, g3list, s3list = [], [], [], []
        try:
            dt = self._read_data()
        except Exception as e:
            raise command.ExecutionError("Import failure", self, e)
        if "c2" in dt:
            for i in range(len(dt["c2"])):
                c2list.append((self._c2name(i), dt["c2"][i]))
        if "g2" in dt:
            for i in range(len(dt["g2"])):
                g2list.append((self._g2name(i), dt["g2"][i]))
        if "s3" in dt:
            for i in range(len(dt["s3"])):
                s3list.append((self._s3name(i), dt["s3"][i]))
        if "g3" in dt:
            for i in range(len(dt["g3"])):
                g3list.append((self._g3name(i), dt["g3"][i]))
        return g2list, [], c2list, [], g3list, [], s3list, []

    def _g2name(self, i):
        if self.options['name'] == '':
            return "Grid1"
        else:
            return self.options['name']

    def _g3name(self, i):
        if self.options['name'] == '':
            return "Grid3D_1"
        else:
            return self.options['name']

    def _c2name(self, i):
        if self.options['name'] == '':
            return "Contour1"
        else:
            return self.options['name']

    def _s3name(self, i):
        if self.options['name'] == '':
            return "Surface1"
        else:
            return self.options['name']

    def __check_file_existance(self):
        'checks self.filename existance and offers to enter new filename'
        fn = self.options['filename']
        if not os.path.isfile(fn):
            raise command.ExecutionError("File not found %s" % fn, self)

    #functions for overriding
    def _read_data(self):
        """ -> {"c2": [list of contours], "g2": [list of grids], "g3":, "s3"}
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
        return {"c2": [cont]}

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
        return {"g2": [grid]}

    #functions for overriding
    def _read_grid(self):
        "-> grid2.Grid2 or raise an Exception"
        raise NotImplementedError


class _AbstractImportGrid3(_AbstractImport):
    'Import contour basic class'
    def __init__(self, argsdict):
        super(_AbstractImportGrid3, self).__init__(argsdict)

    def doc(self):
        return "Import grid3d from %s" % os.path.basename(
            self.options['filename'])

    def _read_data(self):
        #build
        grid = self._read_grid()
        return {"g3": [grid]}

    #functions for overriding
    def _read_grid(self):
        "-> grid3.Grid3 or raise an Exception"
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
        nm, fn, bt, fc = [
            self.options[x] for x in
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


class ImportContoursNative(_AbstractImport):
    'Import all 2d contours from a file'
    def __init__(self, argsdict):
        super(ImportContoursNative, self).__init__(argsdict)

    def doc(self):
        return "Import contours from %s" % os.path.basename(
            self.options['filename'])

    @classmethod
    def _arguments_types(cls):
        """ name - new grids name,
            filename - filename,
            contname - specific name,
            all
        """
        return {'name': command.BasicOption(str),
                'filename': command.BasicOption(str),
                'contname': command.BasicOption(str),
                'all': command.BoolOption()
                }

    def _c2name(self, i):
        return self._c2names[i]

    def _read_data(self):
        so = self.options
        c_reader, c_creader = 0, []
        try:
            c_reader = hmcore.hmxml_read(so['filename'])
            if not so['all']:
                if so['contname'] != '':
                    c_creader = hmcore.hmxml_query(
                        c_reader, 'CONTOUR2D[@name="' + so['contname'] + '"]',
                        required=">0")
                else:
                    c_creader = hmcore.hmxml_query(
                        c_reader, 'CONTOUR2D', required=">0")[:1]
            else:
                c_creader = hmcore.hmxml_query(c_reader, 'CONTOUR2D')

            # construct contours
            conts = [None] * len(c_creader)
            self._c2names = [None] * len(c_creader)
            for i in range(len(c_creader)):
                conts[i], self._c2names[i] =\
                    c2core.contour_from_hmxml(c_reader, c_creader[i])
            return {"c2": conts}
        except Exception:
            raise
        finally:
            for cn in c_creader:
                hmcore.hmxml_free_node(cn)
            hmcore.hmxml_free_node(c_reader) if c_reader != 0 else None


#========================= Import grids
class ImportGridsNative(_AbstractImport):
    'Import all grids from a file'
    def __init__(self, argsdict):
        super(ImportGridsNative, self).__init__(argsdict)

    def doc(self):
        return "Import grids from %s" % os.path.basename(
            self.options['filename'])

    @classmethod
    def _arguments_types(cls):
        """ name - new grids name,
            filename - filename
            gridname
            all
        """
        return {'name': command.BasicOption(str),
                'filename': command.BasicOption(str),
                'gridname': command.BasicOption(str),
                'all': command.BoolOption()
                }

    def _g2name(self, i):
        return self._g2names[i]

    def _read_data(self):
        so = self.options
        c_reader, c_greader = 0, []
        try:
            # find node
            c_reader = hmcore.hmxml_read(so['filename'])
            if not so['all']:
                if so['gridname'] != '':
                    c_greader = hmcore.hmxml_query(
                        c_reader, 'GRID2D[@name="' + so['gridname'] + '"]',
                        required=">0")[:1]
                else:
                    c_greader = hmcore.hmxml_query(
                        c_reader, 'GRID2D', required=">0")[:1]
            else:
                c_greader = hmcore.hmxml_query(c_reader, 'GRID2D')
            # construct grid
            grids = [None] * len(c_greader)
            self._g2names = [None] * len(c_greader)
            for i in range(len(c_greader)):
                grids[i], self._g2names[i] =\
                    g2core.grid_from_hmxml(c_reader, c_greader[i])
            return {"g2": grids}
        except Exception:
            raise
        finally:
            for gr in c_greader:
                hmcore.hmxml_free_node(gr)
            hmcore.hmxml_free_node(c_reader) if c_reader != 0 else None


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
            w1, w2 = d[:k], d[k + 1:]
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
        if 45 not in dtp and 39 not in dtp:
            raise Exception("ZONE section (39, 45) was not found")
        if 45 not in dtp:
            dtp[45] = []
        if 39 in dtp:
            dtp[45] += dtp[39]
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

    def _parse(self):
        """ ->
            {bindex: bname}
            [ ... points ... ],
            [ ... 2d cells ... ],
            [ [p1, p2, btype], ... ]
        """
        r1, r2, r3, r4 = {}, [], [], []
        r2m, r3m = {}, {}
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
                    r2m[int(ds[0]) - 1] = Point2(float(ds[1]), float(ds[2]))
                i += k + 1
            # 3
            if dt[i] == "$Elements":
                k = int(dt[i + 1])
                for d in dt[i + 2:i + 2 + k]:
                    ds = map(int, d.split())
                    index = ds[0] - 1
                    if ds[1] == 1:
                        # edge
                        tp = ds[3] if ds[2] > 0 else 0
                        if tp > 0:
                            r4.append([ds[-2] - 1, ds[-1] - 1, tp])
                    elif ds[1] == 2:
                        # triangle cell
                        r3m[index] = [ds[-3] - 1, ds[-2] - 1, ds[-1] - 1]
                    elif ds[1] == 3:
                        # quad cell
                        r3m[index] = [ds[-4] - 1, ds[-3] - 1,
                                      ds[-2] - 1, ds[-1] - 1]
                    elif ds[1] == 15:
                        # ignore point cell
                        pass
                    else:
                        raise Exception("Not supported gmsh element type")
                i += k + 1
            i += 1

        # point renumbering for successive indexing starting from zero
        need_point_renumbering = len(r2m) != max(r2m.keys()) + 1
        if need_point_renumbering:
            oldnew = {}
            i = 0
            for ind in sorted(r2m.keys()):
                oldnew[ind] = i
                r2.append(r2m[ind])
                i += 1
            for a in r4:
                a[0] = oldnew[a[0]]
                a[1] = oldnew[a[1]]
            for a in r3m.values():
                for i in range(len(a)):
                    a[i] = oldnew[a[i]]
        else:
            for ind in sorted(r2m.keys()):
                r2.append(r2m[ind])

        for ind in sorted(r3m.keys()):
            r3.append(r3m[ind])

        return r1, r2, r3, r4

    def _add_boundaries(self, bt, bedges):
        # collect all boundaries from bedges which not present in the flow
        needed_boundaries = {}
        flowbtypes = self.receiver.get_bnd_types()
        for b in bedges:
            try:
                if b[2] not in needed_boundaries:
                    # if this doesn't raise hence boundary already presents
                    flowbtypes.get(index=b[2])
            except KeyError:
                nm = '-'.join(["gmsh", "boundary", str(b[2])])
                needed_boundaries[b[2]] = nm

        # give them names from bt
        for k in needed_boundaries.keys():
            try:
                needed_boundaries[k] = bt[k]
            except:
                pass

        # add needed boundaries to the flow
        for k in sorted(needed_boundaries.keys()):
            v = needed_boundaries[k]
            com = contcom.EditBoundaryType({"index": k, "name": v})
            self.addcom.append(com)
            self.addcom[-1].do(self.receiver)

    @staticmethod
    def _set_boundaries(grid, bedges):
        # edge index -> boundary index
        imap = {}
        bed1 = grid.boundary_contours()
        bed = []
        for b in bed1:
            bed.extend(b)
        for be in bedges:
            for e in bed:
                edge = grid.edges[e]
                p1 = edge[0]
                p2 = edge[1]
                if (be[0] == p1 and be[1] == p2) or\
                        (be[0] == p2 and be[1] == p1):
                    imap[e] = be[2]
                    break
        # set
        grid.set_edge_bnd(imap)

    #overriding
    def _read_grid(self):
        bt, pts, cls, bedges = self._parse()
        ret = gdata.grid2.Grid2.from_points_cells2(pts, cls)
        self._add_boundaries(bt, bedges)
        self._set_boundaries(ret, bedges)
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


class ImportSurfacesNative(_AbstractImport):
    'Import surfacess from a file'
    def __init__(self, argsdict):
        super(ImportSurfacesNative, self).__init__(argsdict)

    def doc(self):
        return "Import surfaces from %s" % os.path.basename(
            self.options['filename'])

    @classmethod
    def _arguments_types(cls):
        """ name - new surface name,
            filename - filename
            srfname - surface name if not all
            all - import all found contours?
        """
        return {'name': command.BasicOption(str),
                'filename': command.BasicOption(str),
                'srfname': command.BasicOption(str),
                'all': command.BoolOption(),
                }

    def _s3name(self, i):
        return self._s3names[i]

    def _read_data(self):
        so = self.options
        c_reader, c_sreader = 0, []
        try:
            # find node
            c_reader = hmcore.hmxml_read(so['filename'])
            if not so['all']:
                if so['srfname'] != '':
                    c_sreader = hmcore.hmxml_query(
                        c_reader, 'SURFACE3D[@name="' + so['srfname'] + '"]',
                        required=">0")[:1]
                else:
                    c_sreader = hmcore.hmxml_query(
                        c_reader, 'SURFACE3D', required=">0")[:1]
            else:
                c_sreader = hmcore.hmxml_query(c_reader, 'SURFACE3D')
            # construct surface
            srfs = [None] * len(c_sreader)
            self._s3names = [None] * len(c_sreader)
            for i in range(len(c_sreader)):
                srfs[i], self._s3names[i] =\
                    s3core.surface_from_hmxml(c_reader, c_sreader[i])
            return {"s3": srfs}
        except Exception:
            raise
        finally:
            for gr in c_sreader:
                hmcore.hmxml_free_node(gr)
            hmcore.hmxml_free_node(c_reader) if c_reader != 0 else None


class ImportGrids3Native(_AbstractImport):
    'Import all grids from a file'
    def __init__(self, argsdict):
        super(ImportGrids3Native, self).__init__(argsdict)

    def doc(self):
        return "Import 3d grids from %s" % os.path.basename(
            self.options['filename'])

    @classmethod
    def _arguments_types(cls):
        """ name - new grids name,
            filename - filename
            gridname -
            all -
        """
        return {'name': command.BasicOption(str),
                'filename': command.BasicOption(str),
                'gridname': command.BasicOption(str),
                'all': command.BoolOption()
                }

    def _g3name(self, i):
        return self._g3names[i]

    def _read_data(self):
        so = self.options
        c_reader, c_greader = 0, []
        try:
            cb = self.ask_for_callback()
            cb._callback("Loading xml file", "", 0, 0)
            # find node
            c_reader = hmcore.hmxml_read(so['filename'])
            if not so['all']:
                if so['gridname'] != '':
                    c_greader = hmcore.hmxml_query(
                        c_reader, 'GRID3D[@name="' + so['gridname'] + '"]',
                        required=">0")[:1]
                else:
                    c_greader = hmcore.hmxml_query(
                        c_reader, 'GRID3D', required=">0")[:1]
            else:
                c_greader = hmcore.hmxml_query(c_reader, 'GRID3D')
            # construct grid
            grids = [None] * len(c_greader)
            self._g3names = [None] * len(c_greader)
            for i in range(len(c_greader)):
                subcb = cb.subcallback(i, len(c_greader))
                grids[i], self._g3names[i] =\
                    g3core.grid_from_hmxml(c_reader, c_greader[i], subcb)
            return {"g3": grids}
        except Exception:
            raise
        finally:
            for gr in c_greader:
                hmcore.hmxml_free_node(gr)
            hmcore.hmxml_free_node(c_reader) if c_reader != 0 else None


class ImportAllNative(_AbstractImport):
    'Import all data from a file'
    def __init__(self, argsdict):
        super(ImportAllNative, self).__init__(argsdict)

    def doc(self):
        return "Import all hybmesh objects from %s" % os.path.basename(
            self.options['filename'])

    @classmethod
    def _arguments_types(cls):
        """ name - new grids name,
            filename - filename
        """
        return {'name': command.BasicOption(str),
                'filename': command.BasicOption(str),
                }

    def _g2name(self, i):
        return self._gnames[i]

    def _c2name(self, i):
        return self._cnames[i]

    def _g3name(self, i):
        return self._g3names[i]

    def _s3name(self, i):
        return self._s3names[i]

    def _read_data(self):
        import imex
        so = self.options
        c_reader = 0
        try:
            cb = self.ask_for_callback()
            c_reader = hmcore.hmxml_read(so['filename'])
            conts, grids, grids3d, srfs3d,\
                self._cnames, self._gnames, self._g3names, self._s3names =\
                imex.import_all(cb, c_reader)
            return {"c2": conts, "g2": grids, "g3": grids3d, "s3": srfs3d}
        except Exception:
            raise
        finally:
            hmcore.hmxml_free_node(c_reader) if c_reader != 0 else None
