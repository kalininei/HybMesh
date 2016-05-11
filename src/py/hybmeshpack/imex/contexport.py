

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


def tecplot(cont, fname, bt=None):
    out = ['TITLE="2D Contour from HybMesh"',
           'VARIABLES="X" "Y"']
    # --- main contour
    out.extend(['ZONE T="Contour"',
                'ZONETYPE=FELINESEG',
                'N=%i E=%i' % (cont.n_points(), cont.n_edges()),
                'DATAPACKING=BLOCK'
                ])
    for p in cont.points:
        out.append(str(p.x))
    for p in cont.points:
        out.append(str(p.y))
    for e in cont.edges_points():
        out.append('%i %i' % (e[0] + 1, e[1] + 1))
    # --- subcontours
    # assemble
    edges_by_btypes = {}
    for e in cont.edges_points():
        b = e[2]
        if b not in edges_by_btypes:
            edges_by_btypes[b] = [e]
        else:
            edges_by_btypes[b].append(e)
    for b, elist in edges_by_btypes.iteritems():
        if b == 0:
            name = "boundary0"
        else:
            try:
                name = bt.get(index=b).name
            except:
                name = "boundary%i" % b
        out.extend(['ZONE T="%s"' % name,
                    'ZONETYPE=FELINESEG',
                    'E=%i' % len(elist),
                    'D=(1 2)'])
        for e in elist:
            out.append('%i %i' % (e[0] + 1, e[1] + 1))

    _write_list_to_file(out, fname)
