#!/usr/bin/env python
import copy
import math
import itertools
import xml.etree.ElementTree as ET
import hybmeshpack.basic.geom as bgeom
import hybmeshpack.basic.proc as bp
import contour2


class Grid2(bgeom.Point2SetStruct):
    def __init__(self):
        super(Grid2, self).__init__()
        #edges -> points indicies
        self.edges = []
        #cells -> ordered edges indicies
        self.cells = []
        #edge index -> boundary type index
        self.bt = {}
        #contour2.GridContour.
        self.cont = None

    def add_edge_bnd(self, ged):
        """({edge_index: boundary index})
            Adds values to edge boundary prop
        """
        for k, v in ged.items():
            self.bt[k] = v

    def set_edge_bnd(self, ged):
        """({edge_index: boundary index})
            Sets values to edge boundary prop
        """
        self.bt = copy.deepcopy(ged)

    def get_edge_bnd(self, iedge):
        try:
            return self.bt[iedge]
        except KeyError:
            return 0

    def deepcopy(self):
        """ overriden from GeomStruct. returns deepcopied grid """
        #avoiding deepcopy of self.cont
        backup_c = self.cont
        self.cont = None
        ret = super(Grid2, self).deepcopy()
        self.cont = backup_c
        if self.cont is not None:
            contour2.GridContour2.copy_from_source(self, ret)
        return ret

    def reflect(self, p1, p2):
        """ overriden from Point2SetStruct """
        super(Grid2, self).reflect(p1, p2)
        self.cells = [c[::-1] for c in self.cells]

    @staticmethod
    def from_points_cells(npt, ncls, pts, cls):
        """ create grid from row points and cells->nodes
            connectivity given as plain arrays """
        it = iter(pts)
        points = map(bgeom.Point2, it, it)
        cls = bp.multifield_list(cls, "auto")
        return Grid2.from_points_cells2(points, cls)

    @staticmethod
    def from_points_cells2(pts, cls):
        """ create grid from points data pts:
                [bgeom.Point2]
            and cells->nodes connectivity from point indicies cls:
                [[p1 index, p2 index, ...], [p3 index, ....], ...]
            All input data will be deepcopied.
        """
        ret = Grid2()
        ret.points = copy.deepcopy(pts)
        #build cells
        edmap = {}
        for er in cls:
            newcell = []
            for j in range(len(er)):
                nd1, nd2 = er[j - 1], er[j]
                nd = (nd1, nd2) if nd1 < nd2 else (nd2, nd1)
                edind = edmap.get(nd)
                if (edind is None):
                    edind = len(edmap)
                    edmap[nd] = edind
                newcell.append(edind)
            ret.cells.append(newcell)
        #build edges
        ret.edges = [[]] * len(edmap)
        for ed, ind in edmap.items():
            ret.edges[ind] = [ed[0], ed[1]]
        return ret

    @staticmethod
    def from_points_edges(pts, eds):
        """ create grid from points data as
                [bgeom.Point2]
            and point->edge connectivity as
                [[point1, point2, cell_left, cell_right, b_type], ...]
            boundary edges has -1 at cell_left/_right.
        """
        ret = Grid2()
        ret.points = copy.deepcopy(pts)
        #edges -> points
        ret.edges = [[x[0], x[1]] for x in eds]
        #edges -> cells
        raw = []
        for e in eds:
            raw.append(e[2])
            raw.append(e[3])
        #number of cells = maximum cell index + 1
        cells = [[] for i in range(max(raw) + 1)]
        eddir = copy.deepcopy(cells)
        it = iter(raw)
        #write all edges and their directions to cells
        for iedge, (c1, c2) in enumerate(zip(it, it)):
            if c1 >= 0:
                cells[c1].append(iedge)
                eddir[c1].append(True)
            if c2 >= 0:
                cells[c2].append(iedge)
                eddir[c2].append(False)

        #order edges (node by node) and write them to the grid
        ret.cells = [[] for i in range(len(cells))]
        for gcl, rcl in zip(ret.cells, cells):
            gcl.append(rcl.pop(0))
            nextnode = ret.edges[gcl[0]][1]
            while len(rcl) > 0:
                for e in rcl:
                    if nextnode == ret.edges[e][0]:
                        nextnode = ret.edges[e][1]
                        break
                    elif nextnode == ret.edges[e][1]:
                        nextnode = ret.edges[e][0]
                        break
                else:
                    raise Exception("Can not order loaded grid edges")
                gcl.append(e)
                rcl.remove(e)

        #revert edges list for cells with clockwise direction.
        for i in range(len(eddir)):
            if not eddir[i][0]:
                ret.cells[i].reverse()

        ret.build_contour()
        bc2 = ret.boundary_contours()
        bc = []
        for b in bc2:
            bc.extend(b)

        for e in bc:
            if eds[e][4] > 0:
                ret.bt[e] = eds[e][4]

        # #read boundary types
        # for i, e in enumerate(eds):
            # if e[4] > 0:
                # ret.bt[i] = e[4]

        return ret

    def n_edges(self):
        ' -> number of edges '
        return len(self.edges)

    def n_cells(self):
        ' -> number of cells '
        return len(self.cells)

    def cell_types_info(self):
        '->{cell dims: number of such cells}'
        ret = {}
        for k, vit in itertools.groupby(sorted(self.cells, key=len), key=len):
            ret[k] = sum(1 for _ in vit)
        return ret

    def cells_nodes_connect(self):
        """ -> [[n1, n2, n3, n4], [n5, n6, ..., ], ...]

        get nodes indicies for each grid cell
        in a counter-clockwise direction
        """
        ret = [[] for i in range(self.n_cells())]
        for cell, cn in zip(self.cells, ret):
            for i in range(-1, len(cell) - 1):
                e1, e2 = self.edges[cell[i]], self.edges[cell[i + 1]]
                cn.append(e1[0 if e1[0] in e2 else 1])
        return ret

    def edges_cells_connect(self):
        """ -> [[c1, c2], [c3, c4], ...]

        get left and right cell index for each edge.
        if edge is boundary cell index for corresponding side is -1
        """
        #fill returning array with -1
        ret = [[-1, -1] for i in range(self.n_edges())]

        #loop over cells
        for icell, cell in enumerate(self.cells):
            for ne in range(-1, len(cell) - 1):
                #pair of consecutive edges
                e1, e2 = cell[ne], cell[ne + 1]
                #check if inclusion of e1 is positive
                e1_direct = self.edges[e1][1] in self.edges[e2]
                #add to result
                ret[e1][0 if e1_direct else 1] = icell
        return ret

    def central_points(self):
        """ -> [Point2] for each cell """
        ret = []
        for pts in self.cells_nodes_connect():
            ret.append(bgeom.Point2(0, 0))
            for pi in pts:
                p = self.points[pi]
                ret[-1].x += p.x
                ret[-1].y += p.y
            ret[-1].x /= len(pts)
            ret[-1].y /= len(pts)
        return ret

    def build_contour(self):
        'fills self.cont field with contour2.GridContour2 object'
        self.cont = contour2.GridContour2(self)

    def boundary_contours(self):
        """ -> [[e1, e2, e3, ....], [e1, e2, e3, ...], ... ]
            returns indicies of edges which form boundary contours.
            All contours have anti-clockwise direction
        """
        #basic implementation.
        ed_cl = self.edges_cells_connect()
        cand_pos = [(i, self.edges[i])
                    for i, e in enumerate(ed_cl) if e[1] < 0]
        cand_neg = [(i, list(reversed(self.edges[i])))
                    for i, e in enumerate(ed_cl) if e[0] < 0]
        cand = cand_pos + cand_neg
        ret = []
        while len(cand) > 0:
            ed_ind = [cand[0][0]]
            pnt_ind = [cand[0][1][0], cand[0][1][1]]
            del cand[0]
            while pnt_ind[0] != pnt_ind[-1]:
                for i, e in enumerate(cand):
                    if e[1][0] == pnt_ind[-1]:
                        break
                else:
                    raise Exception("Could not detect closed contour")
                ed_ind.append(e[0])
                pnt_ind.append(e[1][1])
                del cand[i]
            ret.append(ed_ind)
        return ret

    # ---- save to xml
    def xml_save(self, xmlnode):
        'save to xml Node'
        #basic info
        xmlnode.attrib["tp"] = self.__class__.__name__
        ET.SubElement(xmlnode, "N_POINTS").text = str(self.n_points())
        ET.SubElement(xmlnode, "N_EDGES").text = str(self.n_edges())
        ET.SubElement(xmlnode, "N_CELLS").text = str(self.n_cells())
        #specific
        self._xml_specific_save(xmlnode)
        #boundary types
        #clean from zeros as its default value
        self.bt = {k: v for k, v in self.bt.items() if v != 0}
        ET.SubElement(xmlnode, "BOUNDARY").text = \
            ' '.join(map(str, bp.dict_to_plain(self.bt)))

    def _xml_specific_save(self, xmlnode):
        'Function for overriding. Is called from xml_save()'
        self._xml_write_nodes(ET.SubElement(xmlnode, "POINTS"))
        self._xml_write_edges(ET.SubElement(xmlnode, "EDGES"))

    def _xml_write_nodes(self, xmlnode):
        ET.SubElement(xmlnode, "COORDS").text = " ".join(map(str, self.points))

    def _xml_write_edges(self, xmlnode):
        #points indicies
        ET.SubElement(xmlnode, "PTS_IND").text = \
            " ".join(map(str, itertools.chain.from_iterable(self.edges)))
        #cells connectivity
        ed_cl = self.edges_cells_connect()
        ET.SubElement(xmlnode, "CONNECT").text = \
            " ".join(map(str, itertools.chain.from_iterable(ed_cl)))

    # ---- load from xml
    @classmethod
    def create_from_xml(cls, xmlnode):
        ' create grid from xml node'
        tp = globals()[xmlnode.attrib['tp']]
        g = tp._specific_create_from_xml(xmlnode)
        bnd = xmlnode.find("BOUNDARY")
        if bnd is not None and bnd.text is not None:
            g.bt = bp.plain_to_dict(map(int, bnd.text.split()))
        g.build_contour()
        return g

    @classmethod
    def _specific_create_from_xml(cls, xmlnode):
        'function for overriding'
        g = cls()
        g._load_nodes(xmlnode.find("POINTS"))
        g._load_edges(xmlnode.find("EDGES"))
        return g

    def _load_nodes(self, xmlnode):
        'fill nodes table from xml node'
        raw_coords = map(float, xmlnode.find("COORDS").text.split())
        it = iter(raw_coords)
        self.points = map(bgeom.Point2, it, it)

    def _load_edges(self, xmlnode):
        'fill edges table from xml node '
        #edges -> points
        raw = map(int, xmlnode.find("PTS_IND").text.split())
        it = iter(raw)
        self.edges = map(list, zip(it, it))
        #edges -> cells
        raw = map(int, xmlnode.find("CONNECT").text.split())
        #number of cells = maximum cell index + 1
        cells = [[] for i in range(max(raw) + 1)]
        eddir = copy.deepcopy(cells)
        it = iter(raw)
        #write all edges and their directions to cells
        for iedge, (c1, c2) in enumerate(zip(it, it)):
            if c1 >= 0:
                cells[c1].append(iedge)
                eddir[c1].append(True)
            if c2 >= 0:
                cells[c2].append(iedge)
                eddir[c2].append(False)

        #order edges (node by node) and write them to the grid
        self.cells = [[] for i in range(len(cells))]
        for gcl, rcl in zip(self.cells, cells):
            gcl.append(rcl.pop(0))
            nextnode = self.edges[gcl[0]][1]
            while len(rcl) > 0:
                for e in rcl:
                    if nextnode == self.edges[e][0]:
                        nextnode = self.edges[e][1]
                        break
                    elif nextnode == self.edges[e][1]:
                        nextnode = self.edges[e][0]
                        break
                else:
                    raise Exception("Can not order loaded grid edges")
                gcl.append(e)
                rcl.remove(e)

        #revert edges list for cells with clockwise direction.
        for i in range(len(eddir)):
            if not eddir[i][0]:
                self.cells[i].reverse()

    def _hash(self):
        # returns float number which is calculated by grid entry
        ret = 0.0
        # edges
        for e in self.edges:
            ret += math.sin(e[0] + 2)
            ret += math.cos(e[1] + 1)
        # cells
        for c in self.cells:
            ret += 14.3 * math.sin(len(c))
            for e in c:
                ret += math.cos(e + 3)
        # points
        for p in self.points:
            ret += 2.0 * math.cos(p.x - 2) - 0.33 * math.sin(p.y + 1)

        return ret


# ============================= Rectanualar structured grids
class TetragonGrid(Grid2):
    """ Structured tetragonal grid

    Let us assume that x,y are the local coordinates of the grid area,
    with the center in bottom left corner of the grid.
    Ordering is as follows:
    - nodes are indexed from bottom left angle in a x-major order
    - horizontal edges are indexed by its left node index
    - vertiacal edges are indexed by its bottom node indes
    - global edges indexing starts with horizontal edges
    - cells are indexed by its bottom left node node index
    """
    def __init__(self, nx, ny):
        super(TetragonGrid, self).__init__()
        self.Nx = nx
        self.Ny = ny
        self.__build_topology()

    def __build_topology(self):
        #horizontal edges
        for j in range(0, self.Ny + 1):
            for i in range(0, self.Nx):
                self.edges.append([self.get_point_index(i, j),
                                   self.get_point_index(i + 1, j)])

        #vertical edges
        for j in range(0, self.Ny):
            for i in range(0, self.Nx + 1):
                self.edges.append([self.get_point_index(i, j),
                                   self.get_point_index(i, j + 1)])

        #cells
        for j in range(0, self.Ny):
            for i in range(0, self.Nx):
                self.cells.append([self.get_hedge_index(i, j),
                                   self.get_vedge_index(i + 1, j),
                                   self.get_hedge_index(i, j + 1),
                                   self.get_vedge_index(i, j)])

    def get_point_index(self, i, j):
        return j * (self.Nx + 1) + i

    def get_hedge_index(self, i, j):
        return j * self.Nx + i

    def get_vedge_index(self, i, j):
        return self.Nx * (self.Ny + 1) + j * (self.Nx + 1) + i

    def get_cell_index(self, i, j):
        return j * self.Nx + i


class UnfRectGrid(TetragonGrid):
    ' Represents uniform rectangular grid with constant step '

    def __init__(self, p0, p1, nx, ny):
        super(UnfRectGrid, self).__init__(nx, ny)
        #points
        hx, hy = (p1.x - p0.x) / nx, (p1.y - p0.y) / ny
        for j in range(0, ny + 1):
            y = p0.y + j * hy
            for i in range(0, nx + 1):
                x = p0.x + i * hx
                self.points.append(bgeom.Point2(x, y))


# ============================= Circular structured grids
class CircularGrid(Grid2):
    """ Topology of the structured circular grid

    Nodes indexing:
        [ (outer arch, angles = [0 -> 2pi]),
          (outer arch - 1, angles = [0 -> 2pi]),
          ....
          (inner arch, angles = [0 -> 2pi]),
          (center node if necessary)]

    Cell index equls its outer node with lower angle (top, right).

    Edge indexing: [ (arch cells), (radius cells) ].
    Arch edge index equals its couter-clockwise node index.
    Raduis edge index equals its outer node index.

    Number of grid primitives depends on the inner cell triangulation option.
    if it is not triangulated:
        number of nodes = Na*Nr
        number of edges = Na*Nr + (Na-1)*Nr
        numner of cells = (Na-1)*Nr+1
    else:
        number of nodes = Na*Nr+1
        number of edges = Na*Nr + Na*Nr
        number of cells = Na*Nr
    """

    def __init__(self, na, nr, is_trian):
        super(CircularGrid, self).__init__()
        self.na = na
        self.nr = nr
        self.is_trian = is_trian
        self.__build_topology()

    def node_index(self, ia, ir):
        ' ia - arch index [0, nr], ir - radius section index [0, na]'
        return ia * self.na + ir

    def arch_edge_index(self, ia, ir):
        return ia * self.na + ir

    def rad_edge_index(self, ia, ir):
        return self.na * self.nr + ia * self.na + ir

    def cell_index(self, ia, ir):
        if not self.is_train and ia == self.nr - 1:
            return self.n_cells() - 1
        else:
            return ia * self.na + ir

    def __build_topology(self):
        ' fills edges and cells arrays'
        #--- edges
        del self.edges[:]
        #arch edges
        for i in range(self.nr):
            for j in range(self.na):
                jnext = (j + 1) % self.na
                self.edges.append([self.node_index(i, j),
                                   self.node_index(i, jnext)])
        #radius edges
        for i in range(self.nr - 1):
            for j in range(self.na):
                self.edges.append([self.node_index(i, j),
                                   self.node_index(i + 1, j)])
        if self.is_trian:
            for j in range(self.na):
                self.edges.append([self.node_index(self.nr - 1, j),
                                   self.na * self.nr])
        #--- cells
        del self.cells[:]
        for i in range(self.nr - 1):
            for j in range(self.na):
                jnext = (j + 1) % self.na
                self.cells.append([self.arch_edge_index(i, j),
                                   self.rad_edge_index(i, jnext),
                                   self.arch_edge_index(i + 1, j),
                                   self.rad_edge_index(i, j)])
        if self.is_trian:
            i = self.nr - 1
            for j in range(self.na):
                jnext = (j + 1) % self.na
                self.cells.append([self.arch_edge_index(i, j),
                                   self.rad_edge_index(i, jnext),
                                   self.rad_edge_index(i, j)])
        else:
            e = [self.arch_edge_index(self.nr - 1, j) for j in range(self.na)]
            self.cells.append(e)


class UnfCircGrid(CircularGrid):
    ' Uniform circular grid '
    def __init__(self, p0, rad, na, nr, coef, is_trian):
        super(UnfCircGrid, self).__init__(na, nr, is_trian)
        self.p_center = copy.deepcopy(p0)
        self.rad = rad
        self.coef = coef
        self.acts = []

        #range segmentation using the refinement coef
        rs = bgeom.div_range(0, self.rad, self.nr, self.coef)
        rs.reverse()

        #build points
        for i in range(self.nr):
            r = rs[i]
            for j in range(self.na):
                a = j * 2 * math.pi / self.na
                p = bgeom.Point2(r * math.cos(a) + self.p_center.x,
                                 r * math.sin(a) + self.p_center.y)
                self.points.append(p)
        if self.is_trian:
            self.points.append(copy.deepcopy(self.p_center))


class UnfRingGrid(CircularGrid):
    ' Uniform circular grid '
    def __init__(self, p0, irad, orad, na, nr, coef):
        super(UnfRingGrid, self).__init__(na, nr + 1, False)
        self.cells = self.cells[:-1]
        self.p_center = copy.deepcopy(p0)
        self.irad = irad
        self.orad = orad
        self.coef = coef
        #angle is used for xml input/output of the grid only
        self.acts = []

        #range segmentation using the refinement coef
        rs = bgeom.div_range(self.irad, self.orad, self.nr, self.coef)
        rs.reverse()

        #build points
        for i in range(self.nr):
            r = rs[i]
            for j in range(self.na):
                a = j * 2 * math.pi / self.na
                p = bgeom.Point2(r * math.cos(a) + self.p_center.x,
                                 r * math.sin(a) + self.p_center.y)
                self.points.append(p)
        if self.is_trian:
            self.points.append(copy.deepcopy(self.p_center))

if __name__ == "__main__":
    pass
