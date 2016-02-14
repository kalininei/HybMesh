#!/usr/bin/env python
"2D contours"
import copy
import xml.etree.ElementTree as ET
import HybMeshPyPack.basic.geom as bgeom
import HybMeshPyPack.basic.proc as bp


class AbstractContour2(bgeom.Point2SetStruct):
    "2D abstract contour"
    def __init__(self):
        super(AbstractContour2, self).__init__()

    def bnd_types(self):
        '->set(int). Boundary indicies which present in contour'
        return set([self.edge_bnd(i) for i in range(self.n_edges())])

    def set_edge_bnd(self, edge_bnd):
        "{edge_index -> bnd index}. Reset edge-> boundary dict"
        raise NotImplementedError

    def add_edge_bnd(self, edge_bnd):
        "{edge_index -> bnd index}. Adds key-value to edge-> boundary dict"
        raise NotImplementedError

    def edge_bnd(self, edge_index):
        "->int. Returns edge boundary index"
        raise NotImplementedError

    def n_edges(self):
        '->int. Number of edges'
        raise NotImplementedError

    def edges_points(self):
        """-> [[i1, i2, bnd1], [], ...]  array of edges for contour:
            index_start_point, index_end_point, boundary_index
        """
        raise NotImplementedError

    def xml_save(self, xmlnode):
        'save to xml Node'
        xmlnode.attrib["tp"] = self.__class__.__name__

    @classmethod
    def create_from_xml(cls, xmlnode):
        'create from xml Node'
        tp = globals()[xmlnode.attrib['tp']]
        return tp.create_from_xml(xmlnode)


class Contour2(AbstractContour2):
    '2d contours with independent geometry'
    def __init__(self):
        super(Contour2, self).__init__()
        self.bnds = {}
        self.edges = []
        self._pnts_ind = {}
        self._sorted_edges = None

    def _geom_changed(self):
        self._sorted_edges = None
        super(Contour2, self)._geom_changed()

    def _append_point(self, p):
        'appends point to the end of list'
        if p not in self.points:
            self._pnts_ind[p] = self.n_points()
            self.points.append(p)

    def append_edges(self, eds, btps=[]):
        """ ([[Point, Point], ... ] eds, [boundary indicies])
            adds edges to contour.
            Points are added as is (no deepcopy)
        """
        #add bnd types
        if len(btps) == len(eds):
            for i, b in enumerate(btps):
                self.bnds[self.n_edges() + i] = b
        #add points
        pt = set()
        map(pt.update, eds)
        map(self._append_point, pt)
        #add edges
        self.edges.extend(eds)
        self._geom_changed()

    def n_edges(self):
        return len(self.edges)

    def set_edge_bnd(self, edge_bnd):
        self.bnds = edge_bnd

    def add_edge_bnd(self, edge_bnd):
        for k, v in edge_bnd.items():
            self.bnds[k] = v

    def edge_bnd(self, edge_index):
        if edge_index in self.bnds:
            return self.bnds[edge_index]
        else:
            return 0

    def edge_info(self, edge_index):
        '->[point1 index, point2 index, btype]'
        e = self.edges[edge_index]
        return [self._pnts_ind[e[0]], self._pnts_ind[e[1]],
                self.edge_bnd(edge_index)]

    def sorted_edges(self):
        """ -> [[e1, e2, e3, ...], []] - indicies of connected edges.
            direction is arbitrary
        """
        if self._sorted_edges is None:
            self._sorted_edges = []
            ep = range(self.n_edges())
            while len(ep) > 0:
                eds = [ep.pop(0)]
                nextp = self.edges[eds[0]][1]
                while 1:
                    ret = bp.find(lambda i: nextp in self.edges[i], ep)
                    if ret is not None:
                        e = self.edges[ret]
                        nextp = e[0 if e[1] == nextp else 1]
                        ep.remove(ret)
                        eds.append(ret)
                    else:
                        break
                self._sorted_edges.append(eds)
        return copy.deepcopy(self._sorted_edges)

    def edges_points(self):
        """-> [[i1, i2, bnd1], [], ...]  array of edges for contour:
            index_start_point, index_end_point, boundary_index
        """
        return map(self.edge_info, range(self.n_edges()))

    def add_from_abstract(self, acont):
        """Add edges from another abstract contour acont"""
        ep = acont.edges_points()
        epe = [[acont.points[e[0]], acont.points[e[1]]] for e in ep]
        epe = copy.deepcopy(epe)
        epb = [e[2] for e in ep]
        self.append_edges(epe, epb)

    def separate(self):
        """ ->[Contour2]
            build set of single connected contours.
            Returns None if self is single connected by itself
        """
        srt_eds = self.sorted_edges()
        if len(srt_eds) == 1:
            return None
        ret = []
        for srt in srt_eds:
            c = Contour2()
            edg = [self.edges[i] for i in srt]
            bnd = [self.edge_bnd(i) for i in srt]
            c.append_edges(copy.deepcopy(edg), bnd)
            ret.append(c)
        return ret

    def simplify(self, angle):
        """ -> Contour2
            Returns contour where all angles between each connected edges
            are grater then angle (degree).

            Edges will not be splitted if they have different boundary types.

            Return None if resulting contours equals self
        """
        if angle == 0:
            angle = 1e-6
        angle = angle / 180.0 * bgeom.math.pi
        epts = self.edges_points()
        ret = Contour2()

        def check_angle(i, ed):
            """ checks angle for first point of i-th edge of ed and
                concatenates this edge with the previous one if necessary.
                Returns True if edge was concatenated.
            """
            if ed[i][2] != ed[i - 1][2]:
                return False
            p1 = self.points[ed[i - 1][0]]
            p2 = self.points[ed[i][0]]
            p3 = self.points[ed[i][1]]
            a = bgeom.angle_3pnt(p1, p2, p3)
            if abs(a - bgeom.math.pi) >= angle:
                return False
            else:
                ed[i - 1][1] = ed[i][1]
                ed.pop(i)
                return True

        for srt in self.sorted_edges():
            #sort points
            septs = [epts[i] for i in srt]
            for i in range(1, len(septs)):
                if septs[i][0] != septs[i - 1][1]:
                    septs[i] = [septs[i][1], septs[i][0], septs[i][2]]
            #build simplified contour
            closed = len(srt) > 2 and septs[0][0] == septs[-1][1]
            i = 0 if closed else 1
            while i < len(septs):
                if not check_angle(i, septs):
                    i += 1
            #add to returning contour
            bnds = [k[2] for k in septs]
            septs = [[self.points[k[0]], self.points[k[1]]] for k in septs]
            ret.append_edges(copy.deepcopy(septs), bnds)

        return ret if ret.n_edges() != self.n_edges() else None

    def xml_save(self, xmlnode):
        super(Contour2, self).xml_save(xmlnode)
        #nodes
        ET.SubElement(xmlnode, "COORDS").text = self.points_to_str()

        #edges
        a = []
        [a.extend(x[0:2]) for x in self.edges_points()]
        ET.SubElement(xmlnode, "EDGES").text = ' '.join(map(str, a))

        #boundaries
        a = map(self.edge_bnd, range(self.n_edges()))
        ET.SubElement(xmlnode, "BOUNDARY").text = ' '.join(map(str, a))

    @classmethod
    def create_from_xml(cls, xmlnode):
        ret = cls()
        ret.append_points_from_str(xmlnode.find("COORDS").text)
        ret._pnts_ind = {p: i for i, p in enumerate(ret.points)}
        it = iter(map(int, xmlnode.find("EDGES").text.split()))
        ret.edges = [[ret.points[i], ret.points[j]] for i, j in zip(it, it)]
        bnode = xmlnode.find("BOUNDARY")
        if bnode is not None and bnode.text is not None:
            bs = map(int, bnode.text.split())
            d = {k: v for k, v in enumerate(bs) if v != 0}
            ret.set_edge_bnd(d)
        return ret

    @classmethod
    def create_from_abstract(cls, acont):
        ret = cls()
        ret.points = copy.deepcopy(acont.points)
        ep = acont.edges_points()
        ret.bnds = {i: b[2] for i, b in enumerate(ep) if b[2] != 0}
        ret.edges = [[ret.points[i], ret.points[j]]
                for [i, j, _] in ep]
        ret._pnts_ind = {p: i for i, p in enumerate(ret.points)}
        return ret

    @classmethod
    def create_from_point_set(cls, pts, edpnt, bnd=[]):
        """ ->Contour2 from set of points and edges->points connectivity.
            pts - [Point2, ...]. All points will be deepcopied from this list
            edpnt - [[p1, p2], [p1, p2], ...]. Points indicies for all edges.
            bnd - boundary types for each edge. if [] -> default boundary
        """
        ret = cls()
        ret.points = copy.deepcopy(pts)
        ret.bnds = {i: b for i, b in enumerate(bnd) if b != 0}
        ret.edges = [[ret.points[i], ret.points[j]]
                for [i, j] in edpnt]
        ret._pnts_ind = {p: i for i, p in enumerate(ret.points)}
        return ret

    def deepcopy(self):
        return self.create_from_abstract(self)


class ClosedContour2(AbstractContour2):
    "2D single closed contour"
    def __init__(self):
        super(ClosedContour2, self).__init__()
        self.bnds = {}

    def append_points(self, pts, bnd=[]):
        "adds a points list along with boundary types"
        oldlen = self.n_points()
        for b in bnd:
            self.bnds[oldlen] = b
            oldlen += 1
        self.points.extend(pts)
        self._geom_changed()

    #overriden from AbstractContour2
    def n_edges(self):
        return self.n_points()

    def set_edge_bnd(self, edge_bnd):
        self.bnds = edge_bnd

    def add_edge_bnd(self, edge_bnd):
        for k, v in edge_bnd.items():
            self.bnds[k] = v

    def edge_bnd(self, edge_index):
        try:
            return self.bnds[edge_index]
        except KeyError:
            return 0

    def edges_points(self):
        ret = [[i, i + 1, self.edge_bnd(i)] for i in range(self.n_edges())]
        ret[-1][1] = 0
        return ret

    def xml_save(self, xmlnode):
        super(ClosedContour2, self).xml_save(xmlnode)
        #nodes
        ET.SubElement(xmlnode, "COORDS").text = self.points_to_str()

        #boundaries
        a = map(self.edge_bnd, range(self.n_edges()))
        ET.SubElement(xmlnode, "BOUNDARY").text = ' '.join(map(str, a))

    @classmethod
    def create_from_xml(cls, xmlnode):
        ret = cls()
        #points
        ret.fill_points_from_str(xmlnode.find('COORDS').text)
        #boundaries
        bnode = xmlnode.find("BOUNDARY")
        if bnode is not None and bnode.text is not None:
            bs = map(int, bnode.text.split())
            d = {k: v for k, v in enumerate(bs) if v != 0}
            ret.set_edge_bnd(d)
        return ret


class ContourFromGrid2(AbstractContour2):
    "Contour from the grid2 edges"
    def __init__(self, g2, eds_list=[]):
        "(Grid2 g2, list-of-edges-indexes eds_list)"
        super(ContourFromGrid2, self).__init__()
        self.g2 = g2
        g2.add_subscriber_change_geom(self._geom_changed)
        if len(eds_list) == 0:
            return
        self.grid_edges = eds_list
        eds_pnt = [g2.edges[i] for i in eds_list]
        #fill points array
        if eds_pnt[0][1] in eds_pnt[1]:
            self.ipoints = [eds_pnt[0][0], eds_pnt[0][1]]
        else:
            self.ipoints = [eds_pnt[0][1], eds_pnt[0][0]]
        for e in eds_pnt[1:-1]:
            if e[0] == self.ipoints[-1]:
                self.ipoints.append(e[1])
            else:
                self.ipoints.append(e[0])
        self.points = map(g2.points.__getitem__, self.ipoints)

    def _ed_i2gi(self, ed_ind):
        "->int. Conversation from contour edge index to grid edge index"
        return self.grid_edges[ed_ind]

    @classmethod
    def copy_from_source(cls, src_grid, tar_grid, cont):
        """ ->ContourFromGrid2 object
            Make a copy of grid contour from src_grid for conform tar_grid.
        """
        ret = cls(tar_grid)
        ret.grid_edges = copy.deepcopy(cont.grid_edges)
        ret.ipoints = copy.deepcopy(cont.ipoints)
        ret.points = map(tar_grid.points.__getitem__, ret.ipoints)
        return ret

    #--- overriden from AbstractContour
    def set_edge_bnd(self, edge_bnd):
        "{edge_index -> bnd index}. Reset edge-> boundary dict"
        eb = {self._ed_i2gi(k): v for k, v in edge_bnd}
        self.g2.set_edge_bnd(eb)

    def add_edge_bnd(self, edge_bnd):
        "{edge_index -> bnd index}. Adds key-value to edge-> boundary dict"
        eb = {self._ed_i2gi(k): v for k, v in edge_bnd}
        self.g2.add_edge_bnd(eb)

    def edge_bnd(self, edge_index):
        "->int. Returns edge boundary index"
        try:
            return self.g2.bt[self.grid_edges[edge_index]]
        except KeyError:
            return 0

    def n_edges(self):
        '->int. Number of edges'
        return self.n_points()

    def edges_points(self):
        """-> [[i1, i2, bnd1], [], ...]  array of edges for contour:
            index_start_point, index_end_point, boundary_index
        """
        ret = [[i, i + 1, self.edge_bnd(i)] for i in range(self.n_edges())]
        ret[-1][1] = 0
        return ret


class GridContour2(AbstractContour2):
    "Contour collection from grid2"
    def __init__(self, g2=None):
        super(GridContour2, self).__init__()
        self.g2 = g2
        if g2 is not None:
            self.cnts = [ContourFromGrid2(g2, x)
                    for x in g2.boundary_contours()]
            self.points = sum([x.points for x in self.cnts], [])
            g2.add_subscriber_change_geom(self._geom_changed)

    def _ed_glob2loc(self, ed_ind):
        """ ->(edge_index, ContourFromGrid2).
            Conversation from global contour edge index to
            local subcontour index"""
        for c in self.cnts:
            if ed_ind >= c.n_edges():
                ed_ind -= c.n_edges()
            else:
                return (c, ed_ind)
        else:
            raise KeyError

    def _pts_loc2glob(self, cont_index, pnt_index):
        '->int. Conversation from local point index to self.points index'
        return sum([self.cnts[i].n_points() for i in range(cont_index)]) + \
                pnt_index

    def _ed_i2gi(self, ed_ind):
        "->int. Conversation from contour edge index to grid edge index"
        v = self._ed_glob2loc(ed_ind)
        return v[0]._ed_i2gi(v[1])

    @classmethod
    def copy_from_source(cls, src_grid, tar_grid):
        """ Make a copy of contour of src_grid.cont for tar_grid.
            Fills tar_grid.cont field and returns GridContour2 object
        """
        cont = cls()
        cont.cnts = [ContourFromGrid2.copy_from_source(src_grid, tar_grid, c1)
                for c1 in src_grid.cont.cnts]
        cont.points = sum([x.points for x in cont.cnts], [])
        tar_grid.add_subscriber_change_geom(cont._geom_changed)
        tar_grid.cont = cont
        cont.g2 = tar_grid
        return cont

    #--- overriden from AbstractContour
    def set_edge_bnd(self, edge_bnd):
        "{edge_index -> bnd index}. Reset edge-> boundary dict"
        grid_edge_bnd = {}
        for k, v in edge_bnd.items():
            try:
                grid_edge_bnd[self._ed_i2gi(k)] = v
            except:
                pass
        self.g2.set_edge_bnd(grid_edge_bnd)

    def add_edge_bnd(self, edge_bnd):
        "{edge_index -> bnd index}. Adds key-value to edge-> boundary dict"
        grid_edge_bnd = {}
        for k, v in edge_bnd.items():
            try:
                grid_edge_bnd[self._ed_i2gi(k)] = v
            except:
                pass
        self.g2.add_edge_bnd(grid_edge_bnd)

    def edge_bnd(self, edge_index):
        "->int. Returns edge boundary index"
        try:
            v = self._ed_glob2loc(edge_index)
            return v[0].edge_bnd(v[1])
        except Exception:
            return 0

    def n_edges(self):
        '->int. Number of edges'
        return sum([x.n_edges() for x in self.cnts])

    def edges_points(self):
        """-> [[i1, i2, bnd1], [], ...]  array of edges for contour:
            index_start_point, index_end_point, boundary_index
        """
        ret = []
        for ic, c in enumerate(self.cnts):
            a = c.edges_points()
            for j in range(len(a)):
                a[j][0] = self._pts_loc2glob(ic, a[j][0])
                a[j][1] = self._pts_loc2glob(ic, a[j][1])
            ret.extend(a)
        return ret
