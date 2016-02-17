import xml.etree.ElementTree as ET
import itertools
import writexml


def _write_list_to_file(lst, fname):
    'writes list of strings to file fname'
    lst = [s + '\n' for s in lst]
    with open(fname, 'w') as outf:
        outf.writelines(lst)


def _check_for_34_grid(grid):
    """ raises if grid has non-triangle/tetrahedral cells
        -> (n3, n4). Returns number of 3/4 cells
    """
    cellsinfo = grid.cell_types_info()
    c3 = cellsinfo[3] if 3 in cellsinfo else 0
    c4 = cellsinfo[4] if 4 in cellsinfo else 0
    if c3 + c4 != grid.n_cells():
        raise Exception("Only triangle/tetrahedral grids could be processed")
    return c3, c4


def hmg(grid, fname):
    outp = ET.Element("GRID2")
    grid.xml_save(outp)
    writexml.writexml(outp, fname)


def vtk(grid, fname):
    out = ["# vtk DataFile Version 3.0",
           grid.__class__.__name__, "ASCII"]
    #Points
    out.append("DATASET UNSTRUCTURED_GRID")
    out.append("POINTS %i float" % grid.n_points())
    out.extend(
        ['%s 0' % str(p) for p in grid.points]
    )
    #Cells
    cp = grid.cells_nodes_connect()
    cd = sum(len(c) for c in cp) + len(cp)
    out.append("CELLS %i %i" % (len(cp), cd))
    out.extend(
        ['%i %s' % (len(c), ' '.join(map(str, c))) for c in cp]
    )
    #CellTypes
    out.append("CELL_TYPES %i" % grid.n_cells())
    out.extend(['7'] * grid.n_cells())

    _write_list_to_file(out, fname)


def ggen(grid, fname):
    c3, c4 = _check_for_34_grid(grid)
    #points boundaries
    bnd = [0] * grid.n_points()
    for e in itertools.chain.from_iterable(grid.boundary_contours()):
        bnd[grid.edges[e][0]] = 1
        bnd[grid.edges[e][1]] = 1

    for e, b in grid.bt.items():
        bnd[grid.edges[e][0]] = b + 1
        bnd[grid.edges[e][1]] = b + 1

    out = []
    cn = grid.cells_nodes_connect()

    def write_points(pts):
        out.extend(
            ['%s %i %i' % (str(p), b, i + 1)
             for i, (p, b) in enumerate(zip(pts, bnd))]
        )

    def write_eds3(eds):
        out.extend(
            ['%i %i %i  %i' % (c[0] + 1, c[1] + 1, c[2] + 1, i + 1)
             for i, c in enumerate(eds)]
        )

    if c4 == 0:
        out.append('%i %i %i' % (grid.n_points(), c3, 0))
        write_points(grid.points)
        write_eds3(cn)
    else:
        cells3 = []
        if c3 == 0:
            for c in cn:
                cells3.extend([[c[0], c[1], c[2]], [c[0], c[2], c[3]]])
            out.append('%i %i %i' % (grid.n_points(), 2 * c4, c4))
            write_points(grid.points)
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
            out.append('%i %i %i' % (grid.n_points(), len(cells3), 0))
            write_points(grid.points)
            write_eds3(cells3)

    _write_list_to_file(out, fname)


def msh(grid, fname, btypes=None):
    """ btypes - gdata.btypes.BndTypesList for setting boundary names
    """
    def toh(i):
        '->str. integer to hex format'
        return hex(i)[2:]

    c3, c4 = _check_for_34_grid(grid)
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
          for i, c in enumerate(grid.edges_cells_connect())}
    for k, v in grid.bt.iteritems():
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
        "(10 (0 1 %s 1 2))" % toh(grid.n_points()),
        "(10 (1 1 %s 1 2)(" % toh(grid.n_points()),
    ])
    out.extend(['  %16.12e %16.12e' % (p.x, p.y)
                for p in grid.points])
    out.extend(["))"])
    #faces
    fconnect = grid.edges_cells_connect()
    out.extend([
        """(0 "Faces")""",
        "(13 (0 1 %s 0))" % toh(grid.n_edges()),
    ])
    #interior faces
    out.append("(13 (3 1 %s 2 0)(" % toh(len(tps[0])))
    out.extend(["2 %s %s %s %s" % (
        toh(grid.edges[i][0] + 1),
        toh(grid.edges[i][1] + 1),
        toh(fconnect[i][0] + 1),
        toh(fconnect[i][1] + 1)
    ) for i in tps[0]])
    out.append('))')
    #boundary nodes
    c0 = len(tps[0])
    for i, t in enumerate(tps[1:]):
        c1 = c0 + len(t)
        out.append("(13 (%s %s %s 3 0)(" % (toh(4 + i), toh(c0 + 1), toh(c1)))
        out.extend(["2 %s %s %s %s" % (toh(grid.edges[i][0] + 1),
                                       toh(grid.edges[i][1] + 1),
                                       toh(fconnect[i][0] + 1),
                                       toh(fconnect[i][1] + 1))
                    for i in t])
        out.append('))')
        c0 = c1

    #cells
    out.extend([
        """(0 "Cells")""",
        """(12 (0 1 %s 0))""" % toh(grid.n_cells()),
        """(12 (2 1 %s 1 0)(""" % toh(grid.n_cells()),
    ])
    out.extend([('1' if len(c) == 3 else '3') for c in grid.cells])
    out.append('))')
    #zones
    out.extend([
        """(0 "Zones")""",
        "(45 (2 fluid fluid)())",
        "(45 (3 interior default-interior)())"
    ])

    zonenames = []
    if btypes is not None:
        for ti in tpind[1:]:
            zonenames.append(btypes.get(index=ti).name)
    else:
        for ti in tpind[1:]:
            zonenames.append('wall%i' % ti)

    out.extend([
        "(45 (%s wall %s)())" % (toh(4 + i), zonenames[i])
        for i in range(len(tpind) - 1)]
    )
    _write_list_to_file(out, fname)


def gmsh(grid, fname, btypes=None):
    c3, c4 = _check_for_34_grid(grid)
    out = []
    out.append("$MeshFormat")
    out.append("2.2 0 8")
    out.append("$EndMeshFormat")
    # bc
    out.append("$PhysicalNames")
    bset = set()
    if btypes is not None:
        for v in grid.bt.values():
            bset.add(v)
    out.append(str(len(bset) + 1))
    intphys = max(bset) + 1
    out.append(' '.join(["2", str(intphys), "\"interior\""]))
    for b in bset:
        out.append(' '.join(["1", str(b),
                   '"' + btypes.get(index=b).name + '"']))
    out.append("$EndPhysicalNames")
    # nodes
    out.append("$Nodes")
    out.append(str(grid.n_points()))
    for i, p in enumerate(grid.points):
        out.append(' '.join([str(i + 1), str(p.x), str(p.y), "0.0"]))
    out.append("$EndNodes")
    # cells elements
    out.append("$Elements")
    cn = grid.cells_nodes_connect()
    be1 = grid.boundary_contours()
    be = []
    for b in be1:
        be.extend(b)
    out.append(str(grid.n_cells() + len(be)))
    for i, c in enumerate(grid.cells):
        etp = "2" if len(c) == 3 else "3"
        nind = ' '.join(map(lambda x: str(x + 1), cn[i]))
        out.append(' '.join([str(i + 1), etp, "2", str(intphys),
                   str(intphys), nind]))
    # line elements
    for i, b in enumerate(be):
        n = str(i + 1 + grid.n_cells())
        nind = ' '.join(map(lambda x: str(x + 1), grid.edges[b]))
        out.append(' '.join([n, "1", "2", str(grid.bt[b]),
                   str(grid.bt[b]), nind]))

    out.append("$EndElements")

    _write_list_to_file(out, fname)
