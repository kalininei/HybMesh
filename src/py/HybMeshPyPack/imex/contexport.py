

def _write_list_to_file(lst, fname):
    'writes list of strings to file fname'
    lst = [s + '\n' for s in lst]
    with open(fname, 'w') as outf:
        outf.writelines(lst)


def vtk(cont, fname):
    out = ["# vtk DataFile Version 3.0",
           cont.__class__.__name__, "ASCII"]
    #Points
    out.append("DATASET UNSTRUCTURED_GRID")
    out.append("POINTS %i float" % cont.n_points())
    out.extend(['%s 0' % str(p) for p in cont.points])
    #Edges
    cp = cont.edges_points()
    cd = 3 * len(cp)
    out.append("CELLS %i %i" % (len(cp), cd))
    out.extend(['2 %i %i' % (c[0], c[1]) for c in cp])
    #CellTypes
    out.append("CELL_TYPES %i" % len(cp))
    out.extend(['3'] * len(cp))
    #Boundary types
    out.append("CELL_DATA %i" % len(cp))
    out.append("SCALARS btype int 1")
    out.append("LOOKUP_TABLE default")
    out.extend(map(str, [x[2] for x in cp]))

    _write_list_to_file(out, fname)
