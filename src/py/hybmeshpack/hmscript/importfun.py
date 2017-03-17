"importing functions"
from hybmeshpack import com
from hybmeshpack.imex import flow_import
from hybmeshpack.hmscript import flow, hmscriptfun
from datachecks import (icheck, Bool, String, ExistingFile)


# Importing grids
@hmscriptfun
def import_grid_hmg(fname, gridname="", allgrids=False):
    """Imports grid from native hmg file.

    :param str fname: file name.

    :param str gridname: name of a grid to find in the file.
      If **gridname** is not specified then the first grid will be taken.

    :param bool allgrids: if True then all grids from the file
      will be imported.

    :returns: grid identifier or list of identifiers if **allgrids** is set.

    """
    icheck(0, ExistingFile())
    icheck(1, String())
    icheck(2, Bool())

    c = com.imcom.ImportGridsNative({
        "filename": fname, "gridname": gridname, "all": allgrids})
    flow.exec_command(c)
    if not allgrids:
        return c.added_grids2()[0]
    else:
        return c.added_grids2()


@hmscriptfun
def import_grid_msh(fname):
    """Imports grid from fluent msh file.

       :param str fname: file name

       :returns: grid identifier

    """
    icheck(0, ExistingFile())

    c = com.imcom.ImportGridMSH({"filename": fname})
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def import_grid_gmsh(fname):
    """Imports grid from gmsh ascii file.

    :param str fname: file name

    :return: grid identifier

    Only triangle and quad elements are supported.

    Boundary types could be exported by passing boundary edges as
    certain Elements of "Line" type. Their physical entity tag
    (first one amoung real tags) will be treated as their boundary index.
    For each such index the new boundary type will be registered in the program
    flow if it has not been registered yet. Name for the new boundary
    type will be taken from PhysicalNames field if it exists, otherwise
    the default name "gmsh-boundary-index" will be used.
    """
    icheck(0, ExistingFile())

    c = com.imcom.ImportGridGMSH({"filename": fname})
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def import3d_grid_hmg(fname, gridname="", allgrids=False):
    """Imports grid from native hmg file.

    :param str fname: file name.

    :param str gridname: name of the grid to import. If not specified
      then the first grid will be loaded.

    :param bool allgrids: if True then all grids from the file
      will be imported.

    :returns: grid identifier/list of identifiers if **allgrids** is set.
    """
    icheck(0, ExistingFile())
    icheck(1, String())
    icheck(2, Bool())

    c = com.imcom.ImportGrids3Native({
        "filename": fname, "gridname": gridname, "all": allgrids})
    flow.exec_command(c)
    if allgrids:
        return c.added_grids3()
    else:
        return c.added_grids3()[0]


# Importing contours
@hmscriptfun
def import_contour_ascii(fname, wbtype=False, force_closed=False):
    """Imports singly connected contour as a sequence of points
    from ascii file.

    :param str fname: file name

    :param bool wbtype: if true reads (x, y, btype_index) for each node
       else reads only coordinates (x, y)

    :param bool force_closed: treat contour as closed even
       if end points are not equal

    :returns: contour identifier
    """
    icheck(0, ExistingFile())
    icheck(1, Bool())
    icheck(2, Bool())
    c = com.imcom.ImportContourASCII({"filename": fname,
                                      "btypes": wbtype,
                                      "force_closed": force_closed})
    flow.exec_command(c)
    return c.added_contours2()[0]


@hmscriptfun
def import_contour_hmc(fname, contname="", allconts=False):
    """Imports contour from hybmesh native format.

      :param str fname: filename

      :param str contname: name of the contour to import. If **contname**
         is not specified then the first one will be taken.

      :param bool allconts: if True then all contours from the file
         will be imported. **contname** will be ignored.

      :returns: contour identifier/list of identifiers if **allconts** is set.

    """
    icheck(0, ExistingFile())
    icheck(1, String())
    icheck(2, Bool())

    c = com.imcom.ImportContoursNative(
        {"filename": fname, "contname": contname, "all": allconts})
    flow.exec_command(c)
    if allconts:
        return c.added_contours2()
    else:
        return c.added_contours2()[0]


@hmscriptfun
def import3d_surface_hmc(fname, srfname="", allsurfs=False):
    """Imports surface from hybmesh native format.

      :param str fname: filename

      :param str srfname: name of the surface to import. If **srfname**
         is not specified then the first one will be taken.

      :param bool allsurfs: if True then all surfaces from the file
         will be imported. **srfname** will be ignored.

      :returns: surface identifier or
         list of identifiers if **allsurfs** is set.

    """
    icheck(0, ExistingFile())
    icheck(1, String())
    icheck(2, Bool())

    c = com.imcom.ImportSurfacesNative(
        {"filename": fname, "srfname": srfname, "all": allsurfs})
    flow.exec_command(c)
    ret = c.added_surfaces3()
    return ret if allsurfs else ret[0]


# data
@hmscriptfun
def import_all_hmd(fname):
    """ Imports all geometry objects found in a file of native format.

    :param str fname: filename

    :returns: list of loaded objects as
      ``[[grid2d ids], [contour2d ids], [grid3d ids], [surface3d ids]]``

    """
    icheck(0, ExistingFile())

    c = com.imcom.ImportAllNative({"filename": fname})
    flow.exec_command(c)
    return [c.added_grids2(), c.added_contours2(), c.added_grids3(),
            c.added_surfaces3()]


@hmscriptfun
def load_project(fname):
    """Loads command flow and data from HybMesh project file.

    :param str fname: file name

    All existing data will be lost.

    See :ref:`hmp-file` for description.
    """
    icheck(0, ExistingFile())
    cb = flow.interface.ask_for_callback()
    flow.to_zero_state()
    flow_import.flow_and_framework_fromfile(fname, flow, cb)
