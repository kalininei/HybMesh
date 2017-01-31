import itertools


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
        raise Exception("Only triangle/quadrangle grids could be processed")
    return c3, c4


def grid2(fname, grid, btypes=None, cb=None):
    c3, c4 = _check_for_34_grid(grid)
    out = []
    out.append("$MeshFormat")
    out.append("2.2 0 8")
    out.append("$EndMeshFormat")
    # needed raw data
    raw_bt = grid.raw_data('btypes')
    raw_vert = grid.raw_data('vertices')
    raw_cellsize = grid.raw_data('cellsizes')
    raw_cellvert = grid.raw_data('cell-vert')
    raw_edgevert = grid.raw_data('edge-vert')
    raw_boundary = grid.raw_data('bedges')
    # extract boundary types
    bdict = {}
    for v in raw_bt:
        if v not in bdict:
            try:
                if v == 0:
                    nm = "default_boundary"
                else:
                    nm = btypes[v]
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
    for i, p in enumerate(raw_vert):
        out.append(' '.join([str(i + 1), str(p[0]), str(p[1]), "0.0"]))
    out.append("$EndNodes")
    # cells elements
    out.append("$Elements")
    it = iter(raw_cellvert)
    for sz in raw_cellsize:
        etp = "2" if sz == 3 else "3"
        sit = itertools.islice(it, sz)
        nind = ' '.join(map(lambda x: str(x + 1), sit))
        out.append(' '.join([str(i + 1), etp, "2", str(intphys),
                   str(intphys), nind]))
    # line elements
    n = grid.n_cells() + 1
    for ib in raw_boundary:
        b, ev = raw_bt[ib], raw_edgevert[ib]
        nind = str(ev[0] + 1) + ' ' + str(ev[1] + 1)
        bind = str(b)
        out.append(' '.join([str(n), "1", "2", bind, bind, nind]))
        n += 1

    out.append("$EndElements")

    _write_list_to_file(out, fname)
