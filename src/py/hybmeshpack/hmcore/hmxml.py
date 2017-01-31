import ctypes as ct
from proc import cport, ccall, ccall_cb, free_cside_array


def open_doc(fname):
    doc, root = ct.c_void_p(0), ct.c_void_p(0)
    ccall(cport.hmxml_open_doc, fname, ct.by_ref(doc), ct.by_ref(root))
    return doc, root


def close_doc(doc, nodes):
    for n in nodes:
        ccall(cport.hmxml_free_node, n)
    ccall(cport.hmxml_free_doc, doc)


def query(node, query_string, required="no"):
    " required = 'no/=0/=1/>0' "
    rnum = ct.c_int(0)
    rans = ct.POINTER(ct.c_void_p)()
    ccall(cport.hmxml_query, node, query, ct.byref(rnum), ct.byref(rans))
    try:
        if required == '>0' and rnum.value == 0:
            raise
        if required == '=1' and rnum.value != 1:
            raise
        if required == '=0' and rnum.value != 0:
            raise
    except:
        for i in range(rnum.value):
            ccall(cport.hmxml_free_node, rans[i])
        raise Exception("Improper number of entries (=%i) for %s" %
                        (rnum.value, query))
    ret = []
    for i in range(rnum.value):
        ret.append(rans[i])
    free_cside_array(rans, 'void*')
    return ret


def doc_to_file(doc, filename):
    if not doc:
        return
    ccall(cport.hmxml_write, doc, filename)


def new_doc():
    rdoc, rroot = ct.c_void_p(), ct.c_void_p()
    ccall(cport.new_writer, ct.byref(rdoc), ct.byref(rroot))


def import_contour2(doc, node, cb=None):
    rcont = ct.c_void_p()
    rname = ct.create_string_buffer(1024)
    ccall(cport.read_contour2, doc, node, ct.byref(rcont), rname)
    return rcont, rname.value


def import_grid2(doc, node, cb=None):
    rgrid = ct.c_void_p()
    rname = ct.create_string_buffer(1024)
    ccall(cport.read_grid2, doc, node, ct.byref(rgrid), rname)
    return rgrid, rname.value


def import_surface3(doc, node, cb=None):
    rsurf = ct.c_void_p()
    rname = ct.create_string_buffer(1024)
    ccall(cport.read_surface3, doc, node, ct.byref(rsurf), rname)
    return rsurf, rname.value


def import_grid3(doc, node, cb=None):
    rgrid = ct.c_void_p()
    rname = ct.create_string_buffer(1024)

    ccall_cb(cb, cport.read_grid3,
             doc, node, ct.byref(rgrid), rname)

    return rgrid, rname.value


def purged_string(doc):
    """ returns string representing hmxml document
        without CONTOUR2D, GRID2D, GRID3D elements
    """
    ret = ct.c_char_p()
    ccall(cport.hmxml_purged_string, ct.byref(ret))
    out = str(ret.value)
    free_cside_array(ret, "char")
    return out
