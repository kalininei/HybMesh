import xml.etree.ElementTree as ET
import itertools
import writexml
from hybmeshpack import hmcore as hmcore
from hybmeshpack.hmcore import g2 as g2core


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


def msh(grid, fname, btypes=None, per_data=None):
    """ btypes - gdata.btypes.BndTypesList for setting boundary names

        per_data: [btype_per, btype_shadow, is_reversed,
                   .....]
    """
    # boundary names
    if btypes is not None:
        c_bnames = hmcore.boundary_names_to_c(btypes)
    else:
        c_bnames = None

    # boundary data
    c_btypes = g2core.boundary_types_to_c(grid)

    # periodic data
    if per_data is not None:
        if len(per_data) % 3 != 0:
            raise Exception("Invalid length of periodic data array")
        periodic_list = []
        for i in range(len(per_data) / 3):
            [tp1, tp2, is_rev] = per_data[3 * i:3 * i + 3]
            if not isinstance(tp1, int) or\
                    not isinstance(tp2, int) or\
                    not isinstance(is_rev, bool):
                raise Exception("Invalid periodic data")
            periodic_list.append(tp1)
            periodic_list.append(tp2)
            periodic_list.append(int(is_rev))
        c_per = hmcore.list_to_c(periodic_list, int)
    else:
        c_per = None

    c_g = g2core.grid_to_c(grid)

    try:
        g2core.to_msh(c_g, fname, c_btypes, c_bnames, c_per)
    finally:
        if c_bnames is not None:
            hmcore.free_boundary_names(c_bnames)
        g2core.free_boundary_types(c_btypes)
        g2core.free_c_grid(c_g)


def gmsh(grid, fname, btypes=None):
    c3, c4 = _check_for_34_grid(grid)
    out = []
    out.append("$MeshFormat")
    out.append("2.2 0 8")
    out.append("$EndMeshFormat")
    # extract boundary types
    bdict = {}
    for v in grid.get_bnd_types():
        if v not in bdict:
            try:
                if v == 0:
                    nm = "default_boundary"
                else:
                    nm = btypes.get(index=v).name
            except:
                nm = "".join(["boundary", str(v)])
            bdict[v] = nm
    intphys = max(bdict.keys()) + 1

    # bc
    out.append("$PhysicalNames")
    out.append(str(len(bdict) + 1))
    out.append(' '.join(["2", str(intphys), "\"interior\""]))
    for b in sorted(bdict.keys()):
        out.append(' '.join(["1", str(b),
                   '"' + bdict[b] + '"']))
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
    n = grid.n_cells() + 1
    for i, b in enumerate(be):
        nind = ' '.join(map(lambda x: str(x + 1), grid.edges[b]))
        bind = str(grid.get_edge_bnd(b))
        out.append(' '.join([str(n), "1", "2", bind, bind, nind]))
        n += 1

    out.append("$EndElements")

    _write_list_to_file(out, fname)


def tecplot(grid, fname, btypes=None):
    c_btypes, c_bnames, c_g = 0, 0, 0
    try:
        c_btypes = g2core.boundary_types_to_c(grid)
        if btypes is not None:
            c_bnames = hmcore.boundary_names_to_c(btypes)
        c_g = g2core.grid_to_c(grid)
        g2core.to_tecplot(c_g, fname, c_btypes, c_bnames)
    except:
        raise
    finally:
        hmcore.free_boundary_names(c_bnames) if c_bnames != 0 else None
        g2core.free_boundary_types(c_btypes) if c_btypes != 0 else None
        g2core.free_c_grid(c_g) if c_g != 0 else None
