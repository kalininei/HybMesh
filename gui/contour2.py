#!/usr/bin/env python
"2D contours"
import copy
import xml.etree.ElementTree as ET
import bp
import bgeom


class AbstractContour2(bgeom.Point2SetStruct):
    "2D abstract contour"
    def __init__(self):
        super(AbstractContour2, self).__init__()

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
    '2d contours whith with independent geometry'
    def __init__(self):
        super(Contour2, self).__init__()
        self.bnds = {}

    def append_points(self, pts, bnd=[]):
        "adds a points list along with boundary types"
        oldlen = self.n_points()
        for b in bnd:
            self.bnds[oldlen] = b
            oldlen += 1
        self.points.extend(pts)
        self._geom_changed()

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

    def xml_save(self, xmlnode):
        super(Contour2, self).xml_save(xmlnode)
        ET.SubElement(xmlnode, "COORDS").text = self.points_to_str()
        a = bp.dict_to_plain(self.bnds)
        ET.SubElement(xmlnode, "BOUNDARY").text = ' '.join(map(str, a))


class ClosedContour2(Contour2):
    "2D contour single closed contour"
    def __init__(self):
        super(ClosedContour2, self).__init__()

    #overriden from AbstractContour2
    def n_edges(self):
        return self.n_points()

    def edges_points(self):
        ret = [[i, i + 1, self.edge_bnd(i)] for i in range(self.n_edges())]
        ret[-1][1] = 0
        return ret

    @classmethod
    def create_from_xml(cls, xmlnode):
        ret = cls()
        #points
        ret.fill_points_from_str(xmlnode.find('COORDS').text)
        #boundaries
        bnd = xmlnode.find('BOUNDARY')
        if bnd is not None:
            bnd = bp.plain_to_dict(map(int, bnd.text.split()))
            ret.set_edge_bnd(bnd)
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
        for k, v in edge_bnd:
            grid_edge_bnd[self._ed_i2gi(k)] = v
        self.g2.set_edge_bnd(grid_edge_bnd)

    def add_edge_bnd(self, edge_bnd):
        "{edge_index -> bnd index}. Adds key-value to edge-> boundary dict"
        grid_edge_bnd = {}
        for k, v in edge_bnd.items():
            grid_edge_bnd[self._ed_i2gi(k)] = v
        self.g2.add_edge_bnd(grid_edge_bnd)

    def edge_bnd(self, edge_index):
        "->int. Returns edge boundary index"
        v = self._ed_glob2loc(edge_index)
        return v[0].edge_bnd(v[1])

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
