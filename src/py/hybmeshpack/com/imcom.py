import os.path
import contcom
from hybmeshpack.hmcore import c2 as c2core
from hybmeshpack.hmcore import hmxml as hmxml
from hybmeshpack.gdata.grid2 import Grid2
from hybmeshpack.gdata.grid3 import Grid3
from hybmeshpack.gdata.srf3 import Surface3
from hybmeshpack.gdata.cont2 import Contour2
import hybmeshpack.imex.native_import as natim
import hybmeshpack.imex.gmsh_import as gmim
import hybmeshpack.imex.fluent_import as flim
import comopt as co
import addremove


# ====================== Abstract Import Commands
class _AbstractImport(addremove.AbstractAddRemove):
    'Import basic class'
    def __init__(self, argsdict):
        super(_AbstractImport, self).__init__(argsdict)

    def _addrem_objects(self):
        self.__check_file_existance()
        g2list, c2list, g3list, s3list = [], [], [], []
        self._init_read()
        c2 = self._read_cont2()
        g2 = self._read_grid2()
        s3 = self._read_surf3()
        g3 = self._read_grid3()
        self._fin_read()
        for i, c in enumerate(c2):
            c2list.append((Contour2(c), self._c2name(i)))
        for i, c in enumerate(g2):
            g2list.append((Grid2(c), self._g2name(i)))
        for i, c in enumerate(s3):
            s3list.append((Surface3(c), self._s3name(i)))
        for i, c in enumerate(g3):
            g3list.append((Grid3(c), self._g3name(i)))

        return g2list, [], c2list, [], g3list, [], s3list, []

    def __check_file_existance(self):
        fn = self.get_option('filename')
        if not os.path.isfile(fn):
            raise Exception("File not found %s" % fn, self)

    #functions for overriding
    def _init_read(self):
        pass

    def _fin_read(self):
        pass

    def _read_grid2(self):
        return []

    def _read_grid3(self):
        return []

    def _read_surf3(self):
        return []

    def _read_cont2(self):
        return []

    def _g2name(self, i):
        return self.get_option('name')

    def _g3name(self, i):
        return self.get_option('name')

    def _c2name(self, i):
        return self.get_option('name')

    def _s3name(self, i):
        return self.get_option('name')


# ========================== Import contours
class ImportContourASCII(_AbstractImport):
    "Import contour from ascii"
    def __init__(self, arg):
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
        return {'name': co.BasicOption(str, None),
                'filename': co.BasicOption(str),
                'simplify': co.BoolOption(False),
                'btypes': co.BoolOption(False),
                'force_closed': co.BoolOption(False)
                }

    def _read_cont2(self):
        fn = self.get_option('filename')
        st = map(float, file(fn).read().split())
        pts, b = [], []
        it = iter(st)
        st = map(float, file(fn).read().split())
        if self.get_option('btypes'):
            for x, y, t in zip(it, it, it):
                pts.append([x, y])
                b.append(int(t))
        else:
            for x, y in zip(it, it):
                pts.append([x, y])
                b.append(0)

        if len(pts) < 2:
            raise Exception("Too few points defined", self)

        return [c2core.build_from_points(
            pts, self.get_option('force_closed'), b)]


class ImportContoursNative(_AbstractImport):
    'Import 2d contours from native file'
    def __init__(self, argsdict):
        super(ImportContoursNative, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        """ name - new grids name,
            filename - filename,
            contname - specific name,
            all
        """
        return {'name': co.BasicOption(str, None),
                'filename': co.BasicOption(str),
                'contname': co.BasicOption(str, ""),
                'all': co.BoolOption(False)
                }

    def _c2name(self, i):
        return self._c2names[i]

    def _read_cont2(self):
        doc, root = 0, 0
        try:
            cb = self.ask_for_callback()
            doc, root = hmxml.open_doc(self.get_option('filename'))
            if self.get_option('all'):
                c2, self._c2names = natim.all_contours2(doc, root, cb)
            else:
                c2, self._c2names = natim.contour2(
                    doc, root, self.get_option('contname'), cb)
                c2, self._c2names = [c2], [self._c2names]
            return c2
        except Exception:
            raise
        finally:
            hmxml.close_doc(doc, [root])


#========================= Import grids
class ImportGridsNative(_AbstractImport):
    'Import grids from native file'
    def __init__(self, argsdict):
        super(ImportGridsNative, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        """ name - new grids name,
            filename - filename
            gridname
            all
        """
        return {'name': co.BasicOption(str, None),
                'filename': co.BasicOption(str),
                'gridname': co.BasicOption(str, ""),
                'all': co.BoolOption(False)
                }

    def _g2name(self, i):
        return self._g2names[i]

    def _read_grid2(self):
        doc, root = 0, 0
        try:
            cb = self.ask_for_callback()
            doc, root = hmxml.open_doc(self.get_option('filename'))
            if self.get_option('all'):
                g2, self._g2names = natim.all_grids2(doc, root, cb)
            else:
                g2, self._g2names = natim.grid2(
                    doc, root, self.get_option('gridname'), cb)
                g2, self._g2names = [g2], [self._g2names]
            return g2
        except Exception:
            raise
        finally:
            hmxml.close_doc(doc, [root])


class _ImportGrid2WithBTypes(_AbstractImport):
    "grid with new boundary names importer: base for fluent and gmsh importers"
    def __init__(self, arg):
        super(_ImportGrid2WithBTypes, self).__init__(arg)
        self.addcom = []  # additional commands to add new boundary types

    def _add_boundaries(self, bt):
        needed_boundaries = {}
        flowbtypes = self.receiver.get_zone_types()
        for k, v in bt.iteritems():
            if k not in flowbtypes:
                needed_boundaries[k] = v

        for k in sorted(needed_boundaries.keys()):
            v = needed_boundaries[k]
            com = contcom.EditBoundaryType({"index": k, "name": v})
            self.addcom.append(com)
            self.addcom[-1].do(self.receiver)

    #overriding
    def _read_grid2(self):
        ret, bt = self._parser(self.get_option('filename'))
        self._add_boundaries(bt)
        return [ret]

    # overriding since additional commands present in reading grid procedure
    def _clear(self):
        super(_ImportGrid2WithBTypes, self)._clear()
        self.addcom = []

    def _undo(self):
        super(_ImportGrid2WithBTypes, self)._undo()
        for c in self.addcom:
            c._undo()

    def _redo(self):
        super(_ImportGrid2WithBTypes, self)._redo()
        for c in self.addcom:
            c._redo()

    # function for overriding
    def _parser(self, fname):
        "-> Grid2.grid2.cdata, {bindex: bname}"
        raise NotImplementedError


class ImportGridMSH(_ImportGrid2WithBTypes):
    "from *.msh format"
    def __init__(self, arg):
        super(ImportGridMSH, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'filename': co.BasicOption(str),
                }

    # overriden
    def _parser(self, fname):
        "-> Grid2.grid2.cdata, {bindex: bname}"
        return flim.grid2(fname)


class ImportGridGMSH(_ImportGrid2WithBTypes):
    "Import grid from msh"
    def __init__(self, arg):
        super(ImportGridGMSH, self).__init__(arg)
        # flow commands which were executed during grid reading
        self.addcom = []

    @classmethod
    def _arguments_types(cls):
        """ name - new grid name,
            filename - filename
        """
        return {'name': co.BasicOption(str, None),
                'filename': co.BasicOption(str),
                }

    # overriden
    def _parser(self, fname):
        "-> Grid2.grid2.cdata, {bindex: bname}"
        return gmim.grid2(fname)


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
        return {'name': co.BasicOption(str, None),
                'filename': co.BasicOption(str),
                'srfname': co.BasicOption(str, ""),
                'all': co.BoolOption(False),
                }

    def _s3name(self, i):
        return self._s3names[i]

    def _read_surf3(self):
        doc, root = 0, 0
        try:
            cb = self.ask_for_callback()
            doc, root = hmxml.open_doc(self.get_option('filename'))
            if self.get_option('all'):
                s3, self._s3names = natim.all_surfaces3(doc, root, cb)
            else:
                s3, self._s3names = natim.surface3(
                    doc, root, self.get_option('srfname'), cb)
                s3, self._s3names = [s3], [self._s3names]
            return s3
        except Exception:
            raise
        finally:
            hmxml.close_doc(doc, [root])


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
        return {'name': co.BasicOption(str, None),
                'filename': co.BasicOption(str),
                'gridname': co.BasicOption(str, ""),
                'all': co.BoolOption(False)
                }

    def _g3name(self, i):
        return self._g3names[i]

    def _read_grid3(self):
        doc, root = 0, 0
        try:
            cb = self.ask_for_callback()
            doc, root = hmxml.open_doc(self.get_option('filename'))
            if self.get_option('all'):
                g3, self._g3names = natim.all_grids3(doc, root, cb)
            else:
                g3, self._g3names = natim.grid3(
                    doc, root, self.get_option('gridname'), cb)
                g3, self._g3names = [g3], [self._g3names]
            return g3
        except Exception:
            raise
        finally:
            hmxml.close_doc(doc, [root])


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
        return {'filename': co.BasicOption(str)}

    def _init_read(self):
        doc, root = 0, 0
        try:
            cb = self.ask_for_callback()
            doc, root = hmxml.open_doc(self.get_option('filename'))
            self.__c2, self.__g2, self.__s3, self.__g3,\
                self._c2names, self._g2names, self._s3names, self._g3names =\
                natim.all_geom(doc, root, cb)
        except Exception:
            raise
        finally:
            hmxml.close_doc(doc, [root])

    def _fin_read(self):
        self.__c2, self.__g2, self.__g3, self.__s3 = None, None, None, None

    def _g2name(self, i):
        return self._g2names[i]

    def _c2name(self, i):
        return self._c2names[i]

    def _g3name(self, i):
        return self._g3names[i]

    def _s3name(self, i):
        return self._s3names[i]

    def _read_grid2(self):
        return self.__g2

    def _read_surf3(self):
        return self.__s3

    def _read_cont2(self):
        return self.__c2

    def _read_grid3(self):
        return self.__g3
