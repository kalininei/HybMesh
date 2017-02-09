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
    raw_btypes = cont.raw_data('btypes')
    raw_vertices = cont.raw_data('vertices')
    raw_edgevert = cont.raw_data('edge-vert')
    # --- main contour
    out.extend(['ZONE T="Contour"',
                'ZONETYPE=FELINESEG',
                'N=%i E=%i' % (cont.n_vertices(), cont.n_edges()),
                'DATAPACKING=BLOCK'
                ])
    for p in raw_vertices:
        out.append(str(p[0]))
    for p in raw_vertices:
        out.append(str(p[1]))
    for e in raw_edgevert:
        out.append('%i %i' % (e[0] + 1, e[1] + 1))
    # --- subcontours
    # assemble
    edges_by_btypes = {}
    for ev, b in zip(raw_edgevert, raw_btypes):
        if b not in edges_by_btypes:
            edges_by_btypes[b] = [list(ev)]
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
