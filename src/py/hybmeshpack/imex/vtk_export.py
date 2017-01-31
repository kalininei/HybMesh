import itertools
from hybmeshpack.hmcore import g3 as g3core


def _write_list_to_file(lst, fname):
    'writes list of strings to file fname'
    lst = [s + '\n' for s in lst]
    with open(fname, 'w') as outf:
        outf.writelines(lst)


def cont2(fname, cont, cb=None):
    out = ["# vtk DataFile Version 3.0",
           cont.__class__.__name__, "ASCII"]
    ne, np = cont.n_edges(), cont.n_points()
    # Raw data
    raw_btypes = cont.raw_data('btypes')
    raw_vert = cont.raw_data('vertices')
    raw_edgevert = cont.raw_data('edge-vert')
    # Points
    out.append("DATASET UNSTRUCTURED_GRID")
    out.append("POINTS %i float" % np)
    out.extend(['%s %s 0' % (str(p[0]), str(p[1])) for p in raw_vert])
    # Edges
    out.append("CELLS %i %i" % (ne, 3 * ne))
    out.extend(['2 %i %i' % (c[0], c[1]) for c in raw_edgevert])
    # CellTypes
    out.append("CELL_TYPES %i" % ne)
    out.extend(['3'] * ne)
    # Boundary types
    out.append("CELL_DATA %i" % ne)
    out.append("SCALARS btype int 1")
    out.append("LOOKUP_TABLE default")
    out.extend(map(str, raw_btypes))

    _write_list_to_file(out, fname)


def grid2(fname, grid, cb=None):
    out = ["# vtk DataFile Version 3.0",
           grid.__class__.__name__, "ASCII"]
    np, nc = grid.n_points(), grid.n_cells()
    # Raw data
    raw_vertices = grid.raw_data('vertices')
    raw_cellsize = grid.raw_data('cellsizes')
    raw_cellvert = grid.raw_data('cell-vert')

    #Points
    out.append("DATASET UNSTRUCTURED_GRID")
    out.append("POINTS %i float" % np)
    out.extend(
        ['%s %s 0' % (str(p[0]), str(p[1])) for p in raw_vertices]
    )
    #Cells
    cd = sum(raw_cellsize) + nc
    out.append("CELLS %i %i" % (nc, cd))
    it = iter(raw_cellvert)
    for sz in raw_cellsize:
        sit = itertools.islice(it, sz)
        out.append('%i %s' % (nc, ' '.join(map(str, sit))))
    #CellTypes
    out.append("CELL_TYPES %i" % nc)
    out.extend(['7'] * nc)

    _write_list_to_file(out, fname)


def grid3(fname, grid, cb=None):
    """ cb -- Callback.CB_CANCEL2 callback object or None
    """
    g3core.to_vtk(grid.cdata, fname, cb)


def vtk_surface(grid, fname, cb=None):
    g3core.surface_to_vtk(grid.cdata, fname, cb)
