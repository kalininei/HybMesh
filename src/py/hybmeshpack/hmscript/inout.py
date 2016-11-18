from hybmeshpack import imex, com
from hybmeshpack.hmscript import flow
from . import ExportError
from . import ExecError


# Exporting grids
def export_grid_vtk(g1, fname):
    """exports grid to vtk format

    Args:
       g1: grid identifier

       fname: output filename
    """
    imex.export_grid("vtk", fname, g1, flow=flow)


def export_grid_hmg(g1, fname, fmt='ascii', afields=[]):
    """exports grid to hybmesh native format

    Args:
      g1: grid identifier

      fname: output filename
    """
    imex.export_grid("hmg", fname, g1, flow=flow,
                     adata={'fmt': fmt, 'afields': afields})


def export_grid_msh(g1, fname, periodic_pairs=None):
    """Exports grid to fluent msh format

    :param g1: 2d grid file identifier

    :param str fname: output filename

    :param list periodic_pairs:
      ``[b-periodic0, b-shadow0, is_reversed0, b-periodic1,
      b-shadow1, is_reversed1, ...]`` list defining periodic boundaries.

      Each periodic condition is defined by three values:

      * ``b-periodic`` - boundary identifier for periodic contour segment
      * ``b-shadow`` - boundary identifier for shadow contour segment
      * ``is_reversed`` - boolean which defines whether shadow contour segment
        should be reversed so that first point of periodic segment be
        equivalent to last point of shadow segment

      Periodic and shadow boundary segments should be singly connected and
      topologically equivalent.

    :raises: hmscript.ExportError

    Only grids with triangle/rectangle cells could be exported.
    """
    try:
        imex.export_grid("msh", fname, g1, flow=flow, adata=periodic_pairs)
    except Exception as e:
        raise ExportError(str(e))


def export_grid_gmsh(g1, fname):
    """exports grid to gmsh ascii format

    :param g1: grid identifier

    :param fname: output filename

    :raises: hmscript.ExportError

    Only grids with triangle/rectangle cells could be exported.

    Boundary edges will be exported as Elements of "Line" type.
    All boundary types which present in grid will be exported as
    Physical Groups with an id
    identical to boundary index and a respective Physical Name.
    """
    try:
        imex.export_grid("gmsh", fname, g1, flow=flow)
    except Exception as e:
        raise ExportError(str(e))


# 3d exports
def export3d_grid_vtk(g1, fname_grid=None, fname_surface=None):
    """Exports 3D grid and its surface to vtk ascii format.

    :param g1: 3D grid file identifier

    :param str-or-None fname_grid: filename for grid output

    :param str-or-None fname_surface: filename for surface output

    :raises: hmscript.ExportError

    Only quadrilateral, wedge and tetrahedron cells could be exported
    as grid. Surface export takes arbitrary grid.

    If a filename is *None* then respective export will be omitted.

    Boundary types are exported as a field called ``boundary_types`` to
    surface output file.
    """
    try:
        if fname_grid is not None:
            imex.export_grid("vtk3d", fname_grid, g1, flow=flow)
        if fname_surface is not None:
            imex.export_grid3_surface("vtk", fname_surface, g1, flow=flow)
    except Exception as e:
        raise ExportError(str(e))


def export3d_grid_msh(g1, fname, periodic_pairs=None):
    """Exports 3D grid to fluent msh ascii format.

    :param g1: 3D grid file identifier

    :param str grid: filename for output

    :param list periodic_pairs:
       ``[periodic-0, shadow-0, periodic-point-0, shadow-point-0,
       periodic-1, shadow-1, periodic-point-1, ...]``

       Each periodic pair is defined by four values:

       * ``periodic`` - boundary identifier for periodic surface
       * ``shadow`` - boundary identifier for shadow surface
       * ``periodic-point`` - point in [x, y, z] format on periodic contour
       * ``shadow-point`` - point in [x, y, z] format on shadow contour

       Given points will be projected to closest vertex on the boundaries
       of respective subsurfaces.

       Periodic and shadow subsurfaces should be singly connected and
       topologically equivalent with respect to given points.
       For surface 2D topology definition periodic/shadow surfaces are taken
       with outside/inside normals respectively.

    :raises: hmscript.ExportError

    """
    try:
        imex.export_grid("msh3d", fname, g1, flow=flow, adata=periodic_pairs)
    except Exception as e:
        raise ExportError(str(e))


def export3d_grid_tecplot(g1, fname):
    """Exports 3D grid to tecplot ascii \*.dat format.

    :param g1: 3D grid file identifier

    :param str grid: filename for output

    :raises: hmscript.ExportError

    A grid zone and zones for each boundary surface defined by boundary type
    will be created in the output file.

    All 3D cells will be saved as FEPOLYHEDRON elements.
    """
    try:
        imex.export_grid("tecplot3d", fname, g1, flow=flow)
    except Exception as e:
        raise ExportError(str(e))


def export3d_grid_hmg(g1, fname, fmt="ascii", afields=[]):
    """ TODO
    """
    try:
        imex.export_grid("hmg3d", fname, g1, flow=flow,
                         adata={'fmt': fmt, 'afields': afields})
    except Exception as e:
        raise ExportError(str(e))


def export_grid_tecplot(g1, fname):
    """exports grid to tecplot ascii \*.dat format

    :param g1: grid identifier

    :param fname: output filename

    :raises: hmscript.ExportError

    All cells will be saved as FEPolygon elements.
    Boundary segments with same boundary type will be converted
    to separate zones.
    """
    try:
        imex.export_grid("tecplot", fname, g1, flow=flow)
    except Exception as e:
        raise ExportError(str(e))


# Exporting contours
def export_contour_vtk(c1, fname):
    """exports contour to vtk format

    :param c1: contour identifier

    :param fname: output filename

    :raises: hmscript.ExportError
    """
    try:
        imex.export_contour("vtk", fname, c1, flow=flow)
    except Exception as e:
        raise ExportError(str(e))


def export_contour_hmc(c1, fname, fmt="ascii"):
    """exports contours to native format

    Args:
       c1: contour identifier

       fname: output filename
    """
    imex.export_contour("hmc", fname, c1, flow=flow, adata={'fmt': fmt})


def export_contour_tecplot(c1, fname):
    """exports contour to tecplot ascii \*.dat format

    :param c1: contour identifier

    :param fname: output filename

    :raises: hmscript.ExportError

    All contour segments will be saved to a zone called "Contour".
    Additional zones will be created for all segments with same
    boundary type.
    """
    try:
        imex.export_contour("tecplot", fname, c1, flow=flow)
    except Exception as e:
        raise ExportError(str(e))


# Importing grids
def import_grid_hmg(fname, gridname=""):
    """Imports grid from native \*.hmg file

    Args:
       fname: file name

    Returns:
       grid identifier
    """
    name = "Grid1"
    if gridname != "":
        name = gridname
    c = com.imcom.ImportGridNative({
        "name": name,
        "filename": fname,
        "gridname": gridname})
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except Exception:
        raise ExecError('import_grid_hmg')


def import_grids_hmg(fname):
    """Imports grid from native \*.hmg file

    Args:
       fname: file name

    Returns:
       grid identifier
    """
    c = com.imcom.ImportGridsNative({"filename": fname})
    try:
        flow.exec_command(c)
        return c._get_added_names()[0]
    except Exception:
        raise ExecError('import_grids_hmg')


def import_grid_msh(fname):
    """Imports grid from fluent \*.msh file

    Args:
       fname: file name

    Returns:
       grid identifier
    """
    c = com.imcom.ImportGridMSH({"filename": fname})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


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
    c = com.imcom.ImportGridGMSH({"filename": fname})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def import3d_grid_hmg(fname, gridname=""):
    """Imports grid from native \*.hmg file

    Args:
       fname: file name

    Returns:
       grid identifier
    """
    name = "Grid1"
    if gridname != "":
        name = gridname
    c = com.imcom.ImportGrid3Native({
        "name": name,
        "filename": fname,
        "gridname": gridname})
    try:
        flow.exec_command(c)
        return c._get_added_names()[2][0]
    except Exception:
        raise ExecError('import3d_grid_hmg')


def import3d_grids_hmg(fname):
    """Imports grid from native \*.hmg file

    Args:
       fname: file name

    Returns:
       grid identifier
    """
    c = com.imcom.ImportGrids3Native({"filename": fname})
    try:
        flow.exec_command(c)
        return c._get_added_names()[2]
    except Exception:
        raise ExecError('import3d_grids_hmg')


# Importing contours
def import_contour_ascii(fname, wbtype=False, force_closed=False):
    """Imports singly connected contour as a sequence of points
    from ascii file.

    Args:
       fname: file name

       wbtype (bool): if true reads (x, y, btype_index) for each node
       else reads only coordinates (x, y)

       force_closed (bool): treat contour as closed even if end points are not
       equal

    Returns:
       contour identifier
    """
    c = com.imcom.ImportContourASCII({"filename": fname,
                                      "btypes": wbtype,
                                      "force_closed": force_closed})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


def import_contour_hmc(fname, contname=""):
    """exports contour to hybmesh native format

    Args:
      g1: grid identifier

      fname: output filename
    """
    c = com.imcom.ImportContourNative(
        {"filename": fname, "contname": contname})
    try:
        flow.exec_command(c)
        return c._get_added_names()[1][0]
    except Exception:
        raise ExecError('import_contour_hmc')


def import_contours_hmc(fname, fmt='ascii'):
    """exports contour to hybmesh native format

    Args:
      g1: grid identifier

      fname: output filename
    """
    c = com.imcom.ImportContoursNative({"filename": fname})
    try:
        flow.exec_command(c)
        return c._get_added_names()[1]
    except Exception:
        raise ExecError('import_contours_hmc')


# data
def import_all_hmd(fname):
    """TODO"""
    c = com.imcom.ImportAllNative({"filename": fname})
    try:
        flow.exec_command(c)
        return c._get_added_names()
    except Exception:
        raise ExecError("import hmd file")


def export_all_hmd(fname, fmt="ascii"):
    """exports contours to native format

    Args:
       c1: contour identifier

       fname: output filename
    """
    try:
        imex.export_all(fname, fmt, flow)
    except Exception as e:
        raise ExportError(str(e))


def save_project(fname):
    """saves current command flow and data to
    HybMesh project file

    Args:
       fname: file name
    """
    imex.write_flow_and_framework_to_file(flow, fname)
