from HybMeshPyPack import imex, com
import HybMeshPyPack.com.imcom
from HybMeshPyPack.hmscript import flow, data


# Exporting grids
def ExportVTK(g1, fname):
    """ exports grid to vtk format
        g1 - grid identifier
        fname - output filename
    """
    g = data.get_grid(name=g1)[2]
    imex.gridexport.vtk(g, fname)


def ExportHMG(g1, fname):
    """ exports grid to hybmesh native format
        g1 - grid identifier
        fname - output filename
    """
    imex.export_grid("hmg", fname, g1, flow=flow)


def ExportMSH(g1, fname):
    """ exports grid to fluent msh format
        g1 - grid identifier
        fname - output filename
    """
    imex.export_grid("msh", fname, g1, flow=flow)


def ExportGMSH(g1, fname):
    """ exports grid to gmsh ascii format
        g1 - grid identifier
        fname - output filename
    """
    imex.export_grid("gmsh", fname, g1, flow=flow)



# Importing grids
def ImportGridHMG(fname):
    """ Imports grid from *.hmg file
        fname - file name
        Returns grid identifier
    """
    c = com.imcom.ImportGridNative({"filename": fname})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def ImportGridMSH(fname):
    """ Imports grid from fluent *.msh file
        fname - file name
        Returns grid identifier
    """
    c = com.imcom.ImportGridMSH({"filename": fname})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def ImportGridGMSH(fname):
    """ Imports grid from gmsh ascii file
        fname - file name
        Returns grid identifier
    """
    c = com.imcom.ImportGridGMSH({"filename": fname})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


# Importing contours
def ImportContourASCII(fname, wbtype=False, force_closed=False):
    """ Imports contour as a sequence of points from ascii file.
        fname - file name
        wbtype - if true reads (x, y, btype_index) for each node
              else reads only coordinates (x, y)
        force_closed - treat contour as closed even if end points are not
                equal
    """
    c = com.imcom.ImportContourASCII({"filename": fname,
        "btypes": wbtype, "force_closed": force_closed})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


# save data
def SaveProject(fname):
    """ saves current command flow and data to
        HybMesh project file
    """
    imex.write_flow_and_framework_to_file(flow, fname)

