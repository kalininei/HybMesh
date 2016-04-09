from hybmeshpack import imex, com
from hybmeshpack.hmscript import flow, data
from . import ExportError


# Exporting grids
def export_grid_vtk(g1, fname):
    """exports grid to vtk format

    Args:
       g1: grid identifier

       fname: output filename
    """
    g = data.get_grid(name=g1)[2]
    imex.gridexport.vtk(g, fname)


def export_grid_hmg(g1, fname):
    """exports grid to hybmesh native format

    Args:
      g1: grid identifier

      fname: output filename
    """
    imex.export_grid("hmg", fname, g1, flow=flow)


def export_grid_msh(g1, fname, periodic_data=None):
    """Exports grid to fluent msh format

    :params g1: 3d grid file identifier

    :params str fname: output filename

    :params periodic_pairs: TODO

    :raises: hmscript.ExportError

    Only grids with triangle/rectangle cells could be exported
    """
    try:
        imex.export_grid("msh", fname, g1, flow=flow, adata=periodic_data)
    except Exception as e:
        raise ExportError(str(e))


def export_grid_gmsh(g1, fname):
    """exports grid to gmsh ascii format

    Args:
       g1: grid identifier

       fname: output filename

    Only grids with triangle/rectangle cells could be exported
    """
    try:
        imex.export_grid("gmsh", fname, g1, flow=flow)
    except Exception as e:
        raise ExportError(str(e))


# 3d exports
def export3d_grid_vtk(g1, fname_grid=None, fname_surface=None):
    """Exports 3d grid to vtk ascii format.

    :params g1: 3d grid file identifier

    :params str-or-None fname_grid: filename for grid output

    :params str-or-None fname_surface: filename for surface output

    If a filename is *None* then respective export will be omitted
    """
    try:
        if fname_grid is not None:
            imex.export_grid("vtk3d", fname_grid, g1, flow=flow)
        if fname_surface is not None:
            imex.export_grid3_surface("vtk", fname_surface, g1, flow=flow)
    except Exception as e:
        raise ExportError(str(e))


def export3d_grid_msh(g1, fname, periodic_pairs=None):
    """Exports 3d grid to fluent msh ascii format.

    :params g1: 3d grid file identifier

    :params str grid: filename for output

    :params periodic_pairs: TODO

    :raises: hmscript.ExportError

    """
    try:
        imex.export_grid("msh3d", fname, g1, flow=flow, adata=periodic_pairs)
    except Exception as e:
        raise ExportError(str(e))


# Exporting contours
def export_contour_vtk(c1, fname):
    """exports contour to vtk format

    Args:
       c1: contour identifier

       fname: output filename
    """
    imex.export_contour("vtk", fname, c1, flow=flow)


# Exporting contours
def export_contour_hmc(c1, fname):
    """exports contour to native format

    Args:
       c1: contour identifier

       fname: output filename
    """
    imex.export_contour("hmc", fname, c1, flow=flow)


# Importing grids
def import_grid_hmg(fname):
    """Imports grid from native \*.hmg file

    Args:
       fname: file name

    Returns:
       grid identifier
    """
    c = com.imcom.ImportGridNative({"filename": fname})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


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
    """Imports grid from gmsh ascii file

    Args:
       fname: file name

    Returns:
       grid identifier
    """
    c = com.imcom.ImportGridGMSH({"filename": fname})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


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


# save data
def save_project(fname):
    """saves current command flow and data to
    HybMesh project file

    Args:
       fname: file name
    """
    imex.write_flow_and_framework_to_file(flow, fname)
