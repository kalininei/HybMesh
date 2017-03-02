from hybmeshpack.hmcore import g2 as g2core
from hybmeshpack.hmcore import g3 as g3core


def _write_list_to_file(lst, fname):
    'writes list of strings to file fname'
    lst = [s + '\n' for s in lst]
    with open(fname, 'w') as outf:
        outf.writelines(lst)


def grid2(fname, grid, btypes={}, cb=None):
    """ btypes: {bindex: bname} """
    g2core.to_tecplot(grid.cdata, fname, btypes, cb)


def grid3(fname, grid, btypes={}, cb=None):
    """ btypes: {bindex: bname} """
    g3core.to_tecplot(grid.cdata, fname, btypes, cb)


def cont2(fname, cont, btypes={}, cb=None):
    out = ['TITLE="2D Contour from HybMesh"',
           'VARIABLES="X" "Y"']
    # --- needed raw data
    raw_btypes = cont.raw_data('bt')
    raw_vertices = cont.raw_data('vert')
    raw_edgevert = cont.raw_data('edge_vert')
    # --- main contour
    out.extend(['ZONE T="Contour"',
                'ZONETYPE=FELINESEG',
                'N=%i E=%i' % (cont.n_vertices(), cont.n_edges()),
                'DATAPACKING=BLOCK'
                ])
    it = iter(raw_vertices)
    for p in zip(it, it):
        out.append(str(p[0]))
    it = iter(raw_vertices)
    for p in zip(it, it):
        out.append(str(p[1]))
    it = iter(raw_edgevert)
    for e in zip(it, it):
        out.append('%i %i' % (e[0] + 1, e[1] + 1))
    # --- subcontours
    # assemble
    edges_by_btypes = {}
    for i, b in enumerate(raw_btypes):
        ev = raw_edgevert[2*i:2*i+2]
        if b not in edges_by_btypes:
            edges_by_btypes[b] = [ev]
        else:
            edges_by_btypes[b].append(e)
    for b, elist in edges_by_btypes.iteritems():
        if b == 0:
            name = "boundary0"
        else:
            try:
                name = btypes[b]
            except KeyError:
                name = "boundary%i" % b
        out.extend(['ZONE T="%s"' % name,
                    'ZONETYPE=FELINESEG',
                    'E=%i' % len(elist),
                    'D=(1 2)'])
        for e in elist:
            out.append('%i %i' % (e[0] + 1, e[1] + 1))
    _write_list_to_file(out, fname)
