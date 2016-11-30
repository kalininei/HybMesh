from hybmeshpack import imex, com
from hybmeshpack.hmscript import flow
from . import ExportError
from . import ExecError


# Exporting grids
def export_grid_vtk(g1, fname):
    """exports grid to vtk format

       :param g1: single or list of grid identifiers

       :param str fname: output filename

       :raises: ExportError
    """
    try:
        imex.export_grid("vtk", fname, g1, flow=flow)
    except Exception as e:
        raise ExportError(str(e))


def export_grid_hmg(g1, fname, fmt='ascii', afields=[]):
    """Exports grid to hybmesh native format.

      :param g1: single or list of grid identifiers

      :param str fname: output filename

      :param str fmt: output data format:

         * ``'ascii'`` - all fields will be saved as text fields,
         * ``'bin'`` - all fields will be saved in binary section,
         * ``'fbin'`` - only floating point fields will be saved
           in binary section.

      :param list-of-str afields: additional data which should be placed
          to output file.

      :raises: hmscript.ExportError

      To save additional data into grid file
      :ref:`user defined fields <udef-fields>`
      place any of these strings into **afields** list:

      * ``'cell-vertices'`` -- cell vertex connectivity. All vertices will
        be written in counterclockwise direction;
      * ``'cell-edges'`` -- cell edge connectivity. All edges will be written
        in counterclockwise direction.

      See :ref:`grid2d-file` for format description.
    """
    try:
        imex.export_grid("hmg", fname, g1, flow=flow,
                         adata={'fmt': fmt, 'afields': afields})
    except Exception as e:
        raise ExportError(str(e))


def export_grid_msh(g1, fname, periodic_pairs=None):
    """Exports grid to fluent msh format

    :param g1: 2d grid file identifier or list of identifiers.

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

    :param g1: single or list of grid identifiers

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


def export_grid_tecplot(g1, fname):
    """exports grid to tecplot ascii \*.dat format

    :param g1: grid identifier or list of identifiers

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


# 3d exports
def export3d_grid_vtk(g1, fname_grid=None, fname_surface=None):
    """Exports 3D grid and its surface to vtk ascii format.

    :param g1: 3D grid file identifier.

    :param str-or-None fname_grid: filename for grid output.

    :param str-or-None fname_surface: filename for surface output.

    :raises: hmscript.ExportError

    Only hexahedron, prism, wedge and tetrahedron cells could be exported
    as a grid. Surface export takes arbitrary grid.

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


def export3d_grid_gmsh(g1, fname):
    """Exports grid to gmsh ascii format.

    :param g1: grid identifier.

    :param str fname: output filename.

    :raises: hmscript.ExportError

    Only grids with tetrahedral/hexahedral/prism/pyramid cells
    could be exported.

    Boundary edges will be exported as Elements of triangle/quadrangle type.
    All boundary types which present in grid will be exported as
    Physical Groups with an id identical to boundary index and
    respective Physical Name.
    """
    try:
        imex.export_grid("gmsh3d", fname, g1, flow=flow)
    except Exception as e:
        raise ExportError(str(e))


def export3d_grid_hmg(g1, fname, fmt="ascii", afields=[]):
    """Exports 3d grid to hybmesh native format.

      :param g1: single or list of grid identifiers

      :param str fname: output filename

      :param str fmt: output data format:

         * ``'ascii'`` - all fields will be saved as text fields,
         * ``'bin'`` - all fields will be saved in binary section,
         * ``'fbin'`` - only floating point fields will be saved
           in binary section.

      :param list-of-str afields: additional data which should be placed
          to output file.

      :raises: hmscript.ExportError

      To save additional data into grid file
      :ref:`user defined fields <udef-fields>`
      place any of these strings into **afields** list:

      * ``'face-vertices'`` - face vertex ordered connectivity,
      * ``'cell-faces'`` - cell face connectivity,
      * ``'cell-vertices'`` - cell vertex connectivity,
      * ``'linfem'`` - tries to write cells-vertex connectivity
        for cell types most widely used in linear fem solvers.
        Supported cell types are: tetrahedron (4 nodes), hexahedron (8),
        prism(6), pyramid(5).
        Record for each of those cells contains points in order prescribed
        by vtk file format. If a cell is not of one of those types
        then a zero length connectivity list will be written for it.

        .. figure:: vtk_cells3d.png
           :width: 500 px

      See :ref:`grid3d-file` for format description.
    """
    try:
        imex.export_grid("hmg3d", fname, g1, flow=flow,
                         adata={'fmt': fmt, 'afields': afields})
    except Exception as e:
        raise ExportError(str(e))


# Exporting contours
def export_contour_vtk(c1, fname):
    """Exports contour to vtk format.

    :param c1: contour identifier or list of identifiers,

    :param str fname: output filename.

    :raises: ExportError
    """
    try:
        imex.export_contour("vtk", fname, c1, flow=flow)
    except Exception as e:
        raise ExportError(str(e))


def export_contour_hmc(c1, fname, fmt="ascii"):
    """Exports contours to native format.

      :param c1: contour identifier or list of identifiers,

      :param str fname: output filename

      :param str fmt: output data format:

         * ``'ascii'`` - all fields will be saved as text fields,
         * ``'bin'`` - all fields will be saved in binary section,
         * ``'fbin'`` - only floating point fields will be saved
           in binary section.

      :raises: ExportError

      See :ref:`contour2d-file` for format description.
    """
    try:
        imex.export_contour("hmc", fname, c1, flow=flow, adata={'fmt': fmt})
    except Exception as e:
        raise ExportError(str(e))


def export_contour_tecplot(c1, fname):
    """Exports contour to tecplot ascii \*.dat format.

    :param c1: contour identifier or list of identifiers,

    :param str fname: output filename.

    :raises: ExportError

    All contour segments will be saved to a zone called "Contour".
    Additional zones will be created for all segments with same
    boundary type.
    """
    try:
        imex.export_contour("tecplot", fname, c1, flow=flow)
    except Exception as e:
        raise ExportError(str(e))


# Importing grids
def import_grid_hmg(fname, gridname="", allgrids=False):
    """Imports grid from native hmg file.

    :param str fname: file name.

    :param str gridname: name of a grid to find in the file.
      If **gridname** is not specified then the first grid will be taken.

    :param bool allgrids: if True then all grids from the file
      will be imported.

    :returns: grid identifier/list of identifiers if **allgrids** is set.

    :raises: ValueError, ExecError.

    """
    c = com.imcom.ImportGridsNative({
        "filename": fname, "gridname": gridname, "all": allgrids})
    try:
        flow.exec_command(c)
        if not allgrids:
            return c._get_added_names()[0][0]
        else:
            return c._get_added_names()[0]
    except Exception:
        raise ExecError('import_grid_hmg')


def import_grid_msh(fname):
    """Imports grid from fluent msh file

       :param str fname: file name

       :returns: grid identifier

       :raises: ExecError
    """
    c = com.imcom.ImportGridMSH({"filename": fname})
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except Exception as e:
        raise ExecError(str(e))


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


def import3d_grid_hmg(fname, gridname="", allgrids=False):
    """Imports grid from native hmg file.

    :param str fname: file name.

    :param str gridname: name of the grid to import. If not specified
      then the first grid will be loaded.

    :param bool allgrids: if True then all grids from the file
      will be imported.

    :returns: grid identifier/list of identifiers if **allgrids** is set.
    """
    c = com.imcom.ImportGrids3Native({
        "filename": fname,
        "gridname": gridname,
        "all": allgrids})
    try:
        flow.exec_command(c)
        if allgrids:
            return c._get_added_names()[2]
        else:
            return c._get_added_names()[2][0]
    except Exception:
        raise ExecError('import3d_grid_hmg')


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


def import_contour_hmc(fname, contname="", allconts=False):
    """Imports contour from hybmesh native format

      :param g1: grid identifier

      :param str fname: output filename

      :param str contname: name of the contour to import. If **contname**
         is not specified then the first one will be taken.

      :param bool allconts: if True then all contours from the file
         will be imported. **contname** will be ignored.

      :returns: contour identifier/list of identifiers if **allconts** is set.

      :raises: ExecError
    """
    c = com.imcom.ImportContoursNative(
        {"filename": fname, "contname": contname, "all": allconts})
    try:
        flow.exec_command(c)
        if allconts:
            return c._get_added_names()[1]
        else:
            return c._get_added_names()[1][0]
    except Exception:
        raise ExecError('import_contour_hmc')


def import3d_surface_hmc(fname, srfname="", allsrfs=False):
    """TODO"""
    c = com.imcom.ImportSurfacesNative(
        {"filename": fname, "srfname": srfname, "all": allsrfs})
    try:
        flow.exec_command(c)
        ret = c._get_added_names()[3]
        if allsrfs:
            return ret[0]
        else:
            return ret
    except Exception:
        raise ExecError('import_contours_hmc')


# data
def import_all_hmd(fname):
    """ Imports all geometry objects found in a file of native format.

    :param str fname: filename

    :returns: list of loaded objects as
      ``[[grid2d ids], [contour2d ids], [grid3d ids], [surface3d ids]]``

    :raises: ExecError
    """
    c = com.imcom.ImportAllNative({"filename": fname})
    try:
        flow.exec_command(c)
        return c._get_added_names()
    except Exception:
        raise ExecError("import hmd file")


def export_all_hmd(fname, fmt="ascii"):
    """Exports all geometrical data to native format.

       :param str fname: output filename

       :param str fmt: output data format:

         * ``'ascii'`` - all fields will be saved as text fields,
         * ``'bin'`` - all fields will be saved in binary section,
         * ``'fbin'`` - only floating point fields will be saved
           in binary section.

       See :ref:`nativeformat` for description.
    """
    try:
        imex.export_all(fname, fmt, flow)
    except Exception as e:
        raise ExportError(str(e))


def save_project(fname, fmt="ascii"):
    """Saves current command flow and data to HybMesh project file.

       :param str fname: file name

       :param str fmt: output data format:

         * ``'ascii'`` - all fields will be saved as text fields,
         * ``'bin'`` - all fields will be saved in binary section,
         * ``'fbin'`` - only floating point fields will be saved
           in binary section.

       :raises: ExportError

       See :ref:`hmp-file` for description.
    """
    try:
        imex.write_flow_and_framework_to_file(flow, fname, fmt)
    except Exception as e:
        raise ExportError(str(e))


def load_project(fname):
    """Loads command flow and data from HybMesh project file.

       :param str fname: file name

       All existing data will be lost.

       See :ref:`hmp-file` for description.
    """
    try:
        flow.to_zero_state()
        imex.read_flow_and_framework_from_file(fname, flow)
    except:
        raise
