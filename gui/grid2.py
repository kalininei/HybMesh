#!/usr/bin/env python
import copy
import math
import bgeom
import xml.etree.ElementTree as ET
import bp


class Grid2(bgeom.GeomStruct):
    def __init__(self):
        super(Grid2, self).__init__()
        #points coords
        self.points = []
        #edges -> points indicies
        self.edges = []
        #cells -> ordered edges indicies
        self.cells = []
        #geometry change substribers:
        #   should contain grid_geom_event(grid2) method
        self._subscribers_change_geom = set()

    def deepcopy(self):
        __tmp = self._subscribers_change_geom
        self._subscribers_change_geom = set()
        ret = super(Grid2, self).deepcopy()
        self._subscribers_change_geom = __tmp
        return ret

    @staticmethod
    def from_points_cells(npt, ncls, pts, cls):
        """ create grid from row points and cells->nodes
            connectivity given as plain arrays """
        ret = Grid2()
        it = iter(pts)
        ret.points = map(bgeom.Point2, it, it)
        #build cells
        edmap = {}
        for er in bp.multifield_list(cls, "auto"):
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

    def n_points(self):
        ' -> number of points '
        return len(self.points)

    def n_edges(self):
        ' -> number of edges '
        return len(self.edges)

    def n_cells(self):
        ' -> number of cells '
        return len(self.cells)

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
                if e1_direct:
                    ret[e1][0] = icell
                else:
                    ret[e1][1] = icell
        return ret

    # ---- events
    def add_subscriber_change_geom(self, subscr):
        """ Adds an object which will get information
            if grid geometry is changed by invocation
            of subsr.grid_geom_event(self) method
        """
        self._subscribers_change_geom.add(subscr)

    def remove_subsriber_change_geom(self, subscr):
        """ Removes subscription from grid geometry event """
        self._subscribers_change_geom.remove(subscr)

    def _geom_changed(self):
        for a in self._subscribers_change_geom:
            a.grid_geom_event(self)

    def bounding_box(self):
        """ -> (bgeom.Point2, bgeom.Point2)

        returns bottom left and top right point of grid bounding box
        """
        x = [p.x for p in self.points]
        y = [p.y for p in self.points]

        return bgeom.Point2(min(x), min(y)), bgeom.Point2(max(x), max(y))

    # ---- Transformations
    def move(self, dx, dy):
        for p in self.points:
            p.x += dx
            p.y += dy
        self._geom_changed()

    def rotate(self, x0, y0, angle):
        self.points = bgeom.rotate_points(self.points, x0, y0, angle)
        self._geom_changed()

    def scale(self, p0, xpc, ypc):
        """ scale the grid using p0 as reference point
            and xpc% and ypc% as scaling procentages
        """
        self.points = bgeom.scale_points(self.points, p0, xpc, ypc)
        self._geom_changed()

    def unscale(self, p0, xpc, ypc):
        """ unscale the grid after scaling with p0 as
            reference point and xpc%, ypc% as scaling procentages
        """
        self.scale(p0, 10000.0 / xpc, 10000.0 / ypc)

    # ---- save to xml
    def xml_save(self, xmlnode):
        ' save to xml Node'
        xmlnode.attrib["tp"] = Factory().string_from_cls(self.__class__)
        ET.SubElement(xmlnode, "N_POINTS").text = str(self.n_points())
        ET.SubElement(xmlnode, "N_EDGES").text = str(self.n_edges())
        ET.SubElement(xmlnode, "N_CELLS").text = str(self.n_cells())
        self._xml_write_nodes(ET.SubElement(xmlnode, "POINTS"))
        self._xml_write_edges(ET.SubElement(xmlnode, "EDGES"))

    def _xml_write_nodes(self, xmlnode):
        ET.SubElement(xmlnode, "COORDS").text = " ".join(map(str, self.points))

    def _xml_write_edges(self, xmlnode):
        #points indicies
        ET.SubElement(xmlnode, "PTS_IND").text = \
                " ".join(map(str, sum(self.edges, [])))
        #cells connectivity
        ed_cl = self.edges_cells_connect()
        ET.SubElement(xmlnode, "CONNECT").text = \
                " ".join(map(str, sum(ed_cl, [])))

    # ---- load from xml
    @classmethod
    def create_from_xml(cls, xmlnode):
        ' create grid from xml node'
        g = cls()
        g._load_nodes(xmlnode.find("POINTS"))
        g._load_edges(xmlnode.find("EDGES"))
        return g

    def _load_nodes(self, xmlnode):
        ' fill nodes table from xml node'
        raw_coords = map(float, xmlnode.find("COORDS").text.split())
        it = iter(raw_coords)
        self.points = map(bgeom.Point2, it, it)

    def _load_edges(self, xmlnode):
        ' fill edges table from xml node '
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

    def cells_nodes_connect(self):
        ret = []
        for j in range(self.Ny):
            for i in range(self.Nx):
                p = self.get_point_index(i, j)
                ret.append([p, p + 1, p + 2 + self.Nx, p + self.Nx + 1])
        return ret

    def get_point_index(self, i, j):
        return j * (self.Nx + 1) + i

    def get_hedge_index(self, i, j):
        return j * self.Nx + i

    def get_vedge_index(self, i, j):
        return self.Nx * (self.Ny + 1) + j * (self.Nx + 1) + i

    def get_cell_index(self, i, j):
        return j * self.Nx + i

    #write to xml
    def _xml_write_nodes(self, xmlnode):
        ET.SubElement(xmlnode, "DIM").text = str(self.Nx) + " " + str(self.Ny)
        super(TetragonGrid, self)._xml_write_nodes(xmlnode)

    def _xml_write_edges(self, xmlnode):
        #no edges is needed since connectivity is postulated
        pass

    #read from xml
    @classmethod
    def create_from_xml(cls, xmlnode):
        [Nx, Ny] = map(int, xmlnode.find("POINTS/DIM").text.split())
        g = cls(Nx, Ny)
        g._loadNodes(xmlnode.find("POINTS"))
        return g


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

        self.angle = 0.0

    def rotate(self, x0, y0, angle):
        self.angle += angle
        super(UnfRectGrid, self).rotate(x0, y0, angle)

    def _xml_write_nodes(self, xmlnode):
        p0 = self.points[self.get_point_index(0, 0)]
        p1 = self.points[self.get_point_index(self.Nx, 0)]
        p2 = self.points[self.get_point_index(0, self.Ny)]
        lx, ly = p0.dist(p1), p0.dist(p2)
        p1 = bgeom.Point2(p0.x + lx, p0.y + ly)
        ET.SubElement(xmlnode, "P0").text = str(p0)
        ET.SubElement(xmlnode, "P1").text = str(p1)
        ET.SubElement(xmlnode, "DIM").text = str(self.Nx) + " " + str(self.Ny)
        ET.SubElement(xmlnode, "ANGLE").text = str(self.angle)

    #read from xml
    @classmethod
    def create_from_xml(cls, xmlnode):
        [Nx, Ny] = map(int, xmlnode.find("POINTS/DIM").text.split())
        p0 = map(float, xmlnode.find("POINTS/P0").text.split())
        p0 = bgeom.Point2(p0[0], p0[1])
        p1 = map(float, xmlnode.find("POINTS/P1").text.split())
        p1 = bgeom.Point2(p1[0], p1[1])
        ret = cls(p0, p1, Nx, Ny)
        an = xmlnode.find("POINTS/ANGLE")
        if an is not None:
            ret.rotate(p0.x, p0.y, float(an.text))
        return ret


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

    def cells_nodes_connect(self):
        ret = []
        # general cells
        for i in range(self.nr - 1):
            for j in range(self.na):
                jnext = (j + 1) % self.na
                ret.append([self.node_index(i, j),
                    self.node_index(i, jnext),
                    self.node_index(i + 1, jnext),
                    self.node_index(i + 1, j)])
        if self.is_trian:
            # inner triangles
            for j in range(self.na):
                jnext = (j + 1) % self.na
                ret.append([self.node_index(self.nr - 1, j),
                    self.node_index(self.nr - 1, jnext),
                    self.n_points() - 1])
        else:
            # inner polygon
            e = [self.node_index(self.nr - 1, j) for j in range(self.na)]
            ret.append(e)

        return ret

    #write to xml
    def _xml_write_nodes(self, xmlnode):
        ET.SubElement(xmlnode, "DIM").text = str(self.na) + " " + str(self.nr)
        super(CircularGrid, self)._xml_write_nodes(xmlnode)

    def _xml_write_edges(self, xmlnode):
        pass

    #read from xml
    @classmethod
    def create_from_xml(cls, xmlnode):
        [na, nr] = map(int, xmlnode.find("POINTS/DIM").text.split())
        g = cls(na, nr)
        g._loadNodes(xmlnode.find("POINTS"))
        return g


class _CircActions:
    def __init__(self, aname, val):
        self.aname = aname
        self.val = val

    def do(self, g):
        if self.aname == "ROTATE":
            g.rotate(g.p_center.x, g.p_center.y, self.val)
        if self.aname == "SCALE":
            g.scale(g.p_center, 100.0, 100.0 * self.val)

    def toxml(self, nd):
        ET.SubElement(nd, self.aname).text = str(self.val)

    @staticmethod
    def fromxml(nd):
        return _CircActions(nd.tag, float(nd.text))


class UnfCircGrid(CircularGrid):
    ' Uniform circular grid '
    def __init__(self, p0, rad, na, nr, coef, is_trian):
        super(UnfCircGrid, self).__init__(na, nr, is_trian)
        self.p_center = p0
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

    def move(self, dx, dy):
        self.p_center.x += dx
        self.p_center.y += dy
        super(UnfCircGrid, self).move(dx, dy)

    def rotate(self, x0, y0, angle):
        [self.p_center] = bgeom.rotate_points([self.p_center], x0, y0, angle)

        if len(self.acts) == 0 or self.acts[-1].aname != "ROTATE":
            self.acts.append(_CircActions("ROTATE", angle))
        else:
            self.acts[-1].val += angle
            if self.acts[-1].val == 0.0:
                self.acts = self.acts[:-1]

        super(UnfCircGrid, self).rotate(x0, y0, angle)

    def scale(self, p0, xpc, ypc):
        [self.p_center] = bgeom.scale_points([self.p_center], p0, xpc, ypc)
        self.rad *= (xpc / 100.0)

        if len(self.acts) == 0 or self.acts[-1].aname != "SCALE":
            self.acts.append(_CircActions("SCALE", ypc / xpc))
        else:
            self.acts[-1].val *= ypc / xpc
            if self.acts[-1].val == 1.0:
                self.acts = self.acts[:-1]

        super(UnfCircGrid, self).scale(p0, xpc, ypc)

    def _xml_write_nodes(self, xmlnode):
        ET.SubElement(xmlnode, "DIM").text = str(self.na) + " " + str(self.nr)
        ET.SubElement(xmlnode, "CENTER").text = str(self.p_center)
        ET.SubElement(xmlnode, "RADIUS").text = str(self.rad)
        ET.SubElement(xmlnode, "RCOEF").text = str(self.coef)
        if len(self.acts) > 0:
            nd = ET.SubElement(xmlnode, "ACTIONS")
            for a in self.acts:
                a.toxml(nd)
        if self.is_trian:
            ET.SubElement(xmlnode, "IN_TRIAN")

    #read from xml
    @classmethod
    def create_from_xml(cls, xmlnode):
        [na, nr] = map(int, xmlnode.find("POINTS/DIM").text.split())
        p0 = map(float, xmlnode.find("POINTS/CENTER").text.split())
        p0 = bgeom.Point2(p0[0], p0[1])
        rad = float(xmlnode.find("POINTS/RADIUS").text)
        coef = float(xmlnode.find("POINTS/RCOEF").text)
        is_trian = xmlnode.find("POINTS/IN_TRIAN") is not None
        ret = cls(p0, rad, na, nr, coef, is_trian)
        anode = xmlnode.findall("POINTS/ACTIONS/*")
        for a in anode:
            act = _CircActions.fromxml(a)
            act.do(ret)
        return ret


class UnfRingGrid(CircularGrid):
    ' Uniform circular grid '
    def __init__(self, p0, irad, orad, na, nr, coef):
        super(UnfRingGrid, self).__init__(na, nr, False)
        self.cells = self.cells[:-1]
        self.p_center = p0
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

    def move(self, dx, dy):
        self.p_center.x += dx
        self.p_center.y += dy
        super(UnfRingGrid, self).move(dx, dy)

    def rotate(self, x0, y0, angle):
        [self.p_center] = bgeom.rotate_points([self.p_center], x0, y0, angle)

        if len(self.acts) == 0 or self.acts[-1].aname != "ROTATE":
            self.acts.append(_CircActions("ROTATE", angle))
        else:
            self.acts[-1].val += angle
            if self.acts[-1].val == 0.0:
                self.acts = self.acts[:-1]

        super(UnfRingGrid, self).rotate(x0, y0, angle)

    def scale(self, p0, xpc, ypc):
        [self.p_center] = bgeom.scale_points([self.p_center], p0, xpc, ypc)
        self.irad *= (xpc / 100.0)
        self.orad *= (xpc / 100.0)

        if len(self.acts) == 0 or self.acts[-1].aname != "SCALE":
            self.acts.append(_CircActions("SCALE", ypc / xpc))
        else:
            self.acts[-1].val *= ypc / xpc
            if self.acts[-1].val == 1.0:
                self.acts = self.acts[:-1]

        super(UnfRingGrid, self).scale(p0, xpc, ypc)

    def _xml_write_nodes(self, xmlnode):
        ET.SubElement(xmlnode, "DIM").text = str(self.na) + " " + str(self.nr)
        ET.SubElement(xmlnode, "CENTER").text = str(self.p_center)
        ET.SubElement(xmlnode, "IRADIUS").text = str(self.irad)
        ET.SubElement(xmlnode, "ORADIUS").text = str(self.orad)
        ET.SubElement(xmlnode, "RCOEF").text = str(self.coef)
        if len(self.acts) > 0:
            nd = ET.SubElement(xmlnode, "ACTIONS")
            for a in self.acts:
                a.toxml(nd)

    #read from xml
    @classmethod
    def create_from_xml(cls, xmlnode):
        [na, nr] = map(int, xmlnode.find("POINTS/DIM").text.split())
        p0 = map(float, xmlnode.find("POINTS/CENTER").text.split())
        p0 = bgeom.Point2(p0[0], p0[1])
        irad = float(xmlnode.find("POINTS/IRADIUS").text)
        orad = float(xmlnode.find("POINTS/ORADIUS").text)
        coef = float(xmlnode.find("POINTS/RCOEF").text)
        ret = cls(p0, irad, orad, na, nr, coef)
        anode = xmlnode.findall("POINTS/ACTIONS/*")
        for a in anode:
            act = _CircActions.fromxml(a)
            act.do(ret)
        return ret


# ============================== grids factory
class Factory(object):
    """
    Represents {grid type -> string code} dictionary
    Factory is used for xml read/write procedures
    """
    instance = None

    def __new__(cls):
        if not cls.instance:
            cls.instance = super(Factory, cls).__new__(cls)
            cls.instance._dict = {
                "Grid2": Grid2,
                "TetragonGrid": TetragonGrid,
                "UnfRectGrid": UnfRectGrid,
                "UnfCircGrid": UnfCircGrid,
                "UnfRingGrid": UnfRingGrid
            }
        return cls.instance

    def string_from_cls(self, cls):
        ind = self._dict.values().index(cls)
        return self._dict.keys()[ind]

    def cls_from_string(self, s):
        return self._dict[s]

    def xml_create(self, xmlnode):
        cls = self.cls_from_string(xmlnode.attrib["tp"])
        return cls.create_from_xml(xmlnode)

if __name__ == "__main__":
    g2 = UnfRectGrid(bgeom.Point2(0, 0), bgeom.Point2(1, 1), 10, 10)
    print g2
    g3 = g2.deepcopy()
    print g3

