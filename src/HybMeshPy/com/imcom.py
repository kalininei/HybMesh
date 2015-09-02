import os.path
import itertools
import basic.geom as bg
import gdata.contour2
import gdata.grid2
import command
import objcom


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
            btypes - read (x y) if false or (x y boundary_type_integer)
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
