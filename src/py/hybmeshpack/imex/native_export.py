from hybmeshpack.hmcore import hmxml
from hybmeshpack.basic.interf import SilentCallbackCancel2
from hybmeshpack.gdata.grid2 import Grid2
from hybmeshpack.gdata.grid3 import Grid3
from hybmeshpack.gdata.cont2 import Contour2
from hybmeshpack.gdata.srf3 import Surface3
from hybmeshpack.hmcore import g2 as g2core
from hybmeshpack.hmcore import g3 as g3core
from hybmeshpack.hmcore import c2 as c2core
from hybmeshpack.hmcore import s3 as s3core


def grid2(doc, node, grid, name, fmt='ascii', afields=[], cb=None):
    g2core.to_hm(doc, node, grid.cdata, name, fmt, afields, cb)


def grid3(doc, node, grid, name, fmt='ascii', afields=[], cb=None):
    g3core.to_hm(doc, node, grid.cdata, name, fmt, afields, cb)


def cont2(doc, node, cont, name, fmt='ascii', cb=None):
    c2core.to_hm(doc, node, cont.cdata, name, fmt, cb)


def surf3(doc, node, surf, name, fmt='ascii', cb=None):
    s3core.to_hm(doc, node, surf.cdata, name, fmt, cb)


def export_all(doc, node, framework, fmt, cb=None):
    if cb is None:
        cb = SilentCallbackCancel2()
    names = framework.get_names()
    for i, nm in enumerate(names):
        cb1 = cb.subcallback(i, len(names))
        ob = framework.get_object(nm)
        if isinstance(ob, Contour2):
            cont2(doc, node, ob, nm, fmt, cb1)
        elif isinstance(ob, Grid2):
            grid2(doc, node, ob, nm, fmt, cb=cb1)
        elif isinstance(ob, Surface3):
            surf3(doc, node, ob, nm, fmt, cb1)
        elif isinstance(ob, Grid3):
            grid3(doc, node, ob, nm, fmt, cb=cb1)


# *_tofile functions
def export_all_tofile(fname, framework, fmt, cb=None):
    doc, root = 0, 0
    try:
        doc, root = hmxml.new_doc()
        export_all(doc, root, framework, fmt, cb)
        hmxml.doc_to_file(doc, fname)
    except:
        raise
    finally:
        hmxml.close_doc(doc, [root])


def grid2_tofile(fname, grids, names, fmt, afields, cb=None):
    if cb is None:
        cb = SilentCallbackCancel2()
    doc, root = 0, 0
    try:
        doc, root = hmxml.new_doc()
        n = min(len(grids), len(names))
        for i in range(n):
            cb1 = cb.subcallback(i, n)
            grid2(doc, root, grids[i], names[i], fmt, afields, cb1)
        hmxml.doc_to_file(doc, fname)
    except:
        raise
    finally:
        hmxml.close_doc(doc, [root])


def grid3_tofile(fname, grids, names, fmt, afields, cb=None):
    if cb is None:
        cb = SilentCallbackCancel2()
    doc, root = 0, 0
    try:
        doc, root = hmxml.new_doc()
        n = min(len(grids), len(names))
        for i in range(n):
            cb1 = cb.subcallback(i, n)
            grid3(doc, root, grids[i], names[i], fmt, afields, cb1)
        hmxml.doc_to_file(doc, fname)
    except:
        raise
    finally:
        hmxml.close_doc(doc, [root])


def cont2_tofile(fname, conts, names, fmt, cb=None):
    if cb is None:
        cb = SilentCallbackCancel2()
    doc, root = 0, 0
    try:
        doc, root = hmxml.new_doc()
        n = min(len(conts), len(names))
        for i in range(n):
            cb1 = cb.subcallback(i, n)
            cont2(doc, root, conts[i], names[i], fmt, cb1)
        hmxml.doc_to_file(doc, fname)
    except:
        raise
    finally:
        hmxml.close_doc(doc, [root])


def surf3_tofile(fname, surfs, names, fmt, cb=None):
    if cb is None:
        cb = SilentCallbackCancel2()
    doc, root = 0, 0
    try:
        doc, root = hmxml.new_doc()
        n = min(len(surfs), len(names))
        for i in range(n):
            cb1 = cb.subcallback(i, n)
            surf3(doc, root, surfs[i], names[i], fmt, cb1)
        hmxml.doc_to_file(doc, fname)
    except:
        raise
    finally:
        hmxml.close_doc(doc, [root])
