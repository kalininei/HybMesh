"file exporting functions"
from hybmeshpack.hmscript import flow, hmscriptfun
from hybmeshpack.imex import vtk_export
from hybmeshpack.imex import native_export
from hybmeshpack.imex import fluent_export
from hybmeshpack.imex import gmsh_export
from hybmeshpack.imex import tecplot_export
from hybmeshpack.imex import flow_export
from datachecks import (icheck, UListOr1, Bool, UList, Point3D,
                        ZType, Grid2D, ACont2D, OneOf, NoneOr,
                        String, Grid3D, ASurf3D, CompoundList)


def _grid2_from_id(gid):
    from hybmeshpack.hmcore import g2
    from hybmeshpack.gdata.grid2 import Grid2
    if not isinstance(gid, list):
        return flow.receiver.get_grid2(gid)
    else:
        ret = map(flow.receiver.get_grid2, gid)
        cd = g2.concatenate([r.cdata for r in ret])
        return Grid2(cd)


def _cont2_from_id(gid):
    from hybmeshpack.hmcore import c2
    from hybmeshpack.gdata.cont2 import Contour2
    if not isinstance(gid, list):
        return flow.receiver.get_any_contour(gid).contour2()
    else:
        ret = map(flow.receiver.get_any_contour, gid)
        ret = map(lambda x: x.contour2(), ret)
        cd = c2.concatenate([r.cdata for r in ret])
        return Contour2(cd)


def _grid3_from_id(gid):
    from hybmeshpack.hmcore import g3
    from hybmeshpack.gdata.grid3 import Grid3
    if not isinstance(gid, list):
        return flow.receiver.get_grid3(gid)
    else:
        ret = map(flow.receiver.get_grid3, gid)
        cd = g3.concatenate([r.cdata for r in ret])
        return Grid3(cd)


def _surf3_from_id(gid):
    from hybmeshpack.hmcore import s3
    from hybmeshpack.gdata.srf3 import Surface3
    if not isinstance(gid, list):
        return flow.receiver.get_any_surface(gid).surface3()
    else:
        ret = map(flow.receiver.get_any_surface, gid)
        ret = map(lambda x: x.surface3(), ret)
        cd = s3.concatenate([r.cdata for r in ret])
        return Surface3(cd)


# Exporting grids
@hmscriptfun
def export_grid_vtk(gid, fname):
    """ Exports 2d grid to vtk format

       :param gid: single or list of 2d grid identifiers

       :param str fname: output filename

       :returns: None
    """
    icheck(0, UListOr1(Grid2D()))
    icheck(1, String())

    grid = _grid2_from_id(gid)
    cb = flow.interface.ask_for_callback()
    vtk_export.grid2(fname, grid, cb)


@hmscriptfun
def export_grid_hmg(gid, fname, fmt='ascii', afields=[]):
    """Exports 2d grid to hybmesh native format.

      :param gid: single or list of grid identifiers

      :param str fname: output filename

      :param str fmt: output data format:

         * ``'ascii'`` - all fields will be saved as text fields,
         * ``'bin'`` - all fields will be saved in binary section,
         * ``'fbin'`` - only floating point fields will be saved
           in binary section.

      :param list-of-str afields: additional data which should be placed
          to output file.

      :returns: None

      To save additional data into grid file
      :ref:`user defined fields <udef-fields>`
      place any of these strings into **afields** list:

      * ``'cell-vertices'`` -- cell vertex connectivity. All vertices will
        be written in counterclockwise direction;
      * ``'cell-edges'`` -- cell edge connectivity. All edges will be written
        in counterclockwise direction.

      See :ref:`grid2d-file` for format description.

    """
    if fmt == "binary":
        fmt = "bin"
    if fmt == "fbinary":
        fmt = "fbin"
    icheck(0, UListOr1(Grid2D()))
    icheck(1, String())
    icheck(2, OneOf('ascii', 'bin', 'fbin'))
    icheck(3, UList(OneOf('cell-vertices', 'cell-edges')))

    names = gid if isinstance(gid, list) else [gid]
    grids = map(flow.receiver.get_grid2, names)
    cb = flow.interface.ask_for_callback()
    native_export.grid2_tofile(fname, grids, names, fmt, afields, cb)


@hmscriptfun
def export_grid_msh(gid, fname, periodic_pairs=[]):
    """Exports grid to fluent msh format

    :param gid: 2d grid file identifier or list of identifiers.

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

    :returns: None

    Only grids with triangle/quadrangle cells could be exported.
    """
    icheck(0, UListOr1(Grid2D()))
    icheck(1, String())
    icheck(2, CompoundList(ZType(), ZType(), Bool()))

    cb = flow.interface.ask_for_callback()
    grid = _grid2_from_id(gid)
    bt = flow.receiver.get_zone_types()
    fluent_export.grid2(fname, grid, bt, periodic_pairs, cb)


@hmscriptfun
def export_grid_gmsh(gid, fname):
    """ Exports grid to gmsh ascii format

    :param gid: single or list of grid identifiers

    :param fname: output filename

    :returns: None

    Only grids with triangle/quadrangle cells could be exported.

    Boundary edges will be exported as Elements of "Line" type.
    All boundary types which present in grid will be exported as
    Physical Groups with an id
    identical to boundary index and a respective Physical Name.
    """
    icheck(0, UListOr1(Grid2D()))
    icheck(1, String())

    cb = flow.interface.ask_for_callback()
    grid = _grid2_from_id(gid)
    bt = flow.receiver.get_zone_types()
    gmsh_export.grid2(fname, grid, bt, cb)


@hmscriptfun
def export_grid_tecplot(gid, fname):
    """exports grid to tecplot ascii \*.dat format

    :param gid: grid identifier or list of identifiers

    :param fname: output filename

    :returns: None

    All cells will be saved as FEPolygon elements.
    Boundary segments with same boundary type will be converted
    to separate zones.
    """
    icheck(0, UListOr1(Grid2D()))
    icheck(1, String())

    cb = flow.interface.ask_for_callback()
    grid = _grid2_from_id(gid)
    bt = flow.receiver.get_zone_types()
    tecplot_export.grid2(fname, grid, bt, cb)


# 3d exports
@hmscriptfun
def export3d_grid_vtk(gid, fname_grid=None, fname_surface=None):
    """Exports 3D grid and its surface to vtk ascii format.

    :param gid: 3D grid file identifier or list of identifiers

    :param str-or-None fname_grid: filename for grid output.

    :param str-or-None fname_surface: filename for surface output.

    Only hexahedron, prism, wedge and tetrahedron cells could be exported
    as a grid. Surface export takes arbitrary grid.

    If a filename is *None* then respective export will be omitted.

    Boundary types are exported as a field called ``boundary_types`` to
    surface output file.
    """
    icheck(0, UListOr1(Grid3D()))
    icheck(1, NoneOr(String()))
    icheck(2, NoneOr(String()))
    if len(fname_grid) == 0:
        fname_grid = None
    if len(fname_surface) == 0:
        fname_surface = None

    cb = flow.interface.ask_for_callback()
    grid = _grid3_from_id(gid)
    if fname_grid is not None:
        vtk_export.grid3(fname_grid, grid, cb)
    if fname_surface is not None:
        vtk_export.grid3_surface(fname_surface, grid, cb)


@hmscriptfun
def export3d_grid_msh(gid, fname, periodic_pairs=[]):
    """Exports 3D grid to fluent msh ascii format.

    :param gid: 3D grid file identifier or list of identifiers

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

    """
    icheck(0, UListOr1(Grid3D()))
    icheck(1, String())
    icheck(2, CompoundList(ZType(), ZType(), Point3D(), Point3D()))

    cb = flow.interface.ask_for_callback()
    grid = _grid3_from_id(gid)
    bt = flow.receiver.get_zone_types()
    fluent_export.grid3(fname, grid, bt, periodic_pairs, cb)


@hmscriptfun
def export3d_grid_tecplot(gid, fname):
    """Exports 3D grid to tecplot ascii \*.dat format.

    :param gid: 3D grid file identifier or list of identifiers

    :param str grid: filename for output

    A grid zone and zones for each boundary surface defined by boundary type
    will be created in the output file.

    All 3D cells will be saved as FEPOLYHEDRON elements.
    """
    icheck(0, UListOr1(Grid3D()))
    icheck(1, String())

    cb = flow.interface.ask_for_callback()
    grid = _grid3_from_id(gid)
    bt = flow.receiver.get_zone_types()
    tecplot_export.grid3(fname, grid, bt, cb)


@hmscriptfun
def export3d_grid_gmsh(gid, fname):
    """Exports 3D grid to gmsh ascii format.

    :param gid: grid identifier or list of identifiers

    :param str fname: output filename.

    Only grids with tetrahedral/hexahedral/prism/pyramid cells
    could be exported.

    Boundary edges will be exported as Elements of triangle/quadrangle type.
    All boundary types which present in grid will be exported as
    Physical Groups with an id identical to boundary index and
    respective Physical Name.
    """
    icheck(0, UListOr1(Grid3D()))
    icheck(1, String())
    cb = flow.interface.ask_for_callback()
    grid = _grid3_from_id(gid)
    bt = flow.receiver.get_zone_types()
    gmsh_export.grid3(fname, grid, bt, cb)


@hmscriptfun
def export3d_grid_hmg(gid, fname, fmt="ascii", afields=[]):
    """Exports 3d grid to hybmesh native format.

      :param gid: single or list of grid identifiers

      :param str fname: output filename

      :param str fmt: output data format:

         * ``'ascii'`` - all fields will be saved as text fields,
         * ``'bin'`` - all fields will be saved in binary section,
         * ``'fbin'`` - only floating point fields will be saved
           in binary section.

      :param list-of-str afields: additional data which should be placed
          to output file.

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
    if fmt == "binary":
        fmt = "bin"
    if fmt == "fbinary":
        fmt = "fbin"
    icheck(0, UListOr1(Grid3D()))
    icheck(1, String())
    icheck(2, OneOf('ascii', 'bin', 'fbin'))
    icheck(3, UList(OneOf('face-vertices', 'cell-faces', 'cell-vertices',
                          'linfem')))
    names = gid if isinstance(gid, list) else [gid]
    grids = map(flow.receiver.get_grid3, names)
    cb = flow.interface.ask_for_callback()
    native_export.grid3_tofile(fname, grids, names, fmt, afields, cb)


# Exporting contours
@hmscriptfun
def export_contour_vtk(cid, fname):
    """Exports contour to vtk format.

    :param cid: contour identifier or list of identifiers,

    :param str fname: output filename.

    :return: None
    """
    icheck(0, UListOr1(ACont2D()))
    icheck(1, String())

    cont = _cont2_from_id(cid)
    cb = flow.interface.ask_for_callback()
    vtk_export.cont2(fname, cont, cb)


@hmscriptfun
def export_contour_hmc(cid, fname, fmt="ascii"):
    """Exports contours to native format.

      :param cid: contour identifier or list of identifiers,

      :param str fname: output filename

      :param str fmt: output data format:

         * ``'ascii'`` - all fields will be saved as text fields,
         * ``'bin'`` - all fields will be saved in binary section,
         * ``'fbin'`` - only floating point fields will be saved
           in binary section.

      :returns: None

      See :ref:`contour2d-file` for format description.
    """
    if fmt == "binary":
        fmt = "bin"
    elif fmt == "fbinary":
        fmt = "fbin"
    icheck(0, UListOr1(ACont2D()))
    icheck(1, String())
    icheck(2, OneOf('ascii', 'bin', 'fbin'))

    names = cid if isinstance(cid, list) else [cid]
    conts = [_cont2_from_id(n) for n in names]
    cb = flow.interface.ask_for_callback()
    native_export.cont2_tofile(fname, conts, names, fmt, cb)


@hmscriptfun
def export_contour_tecplot(cid, fname):
    """Exports contour to tecplot ascii \*.dat format.

    :param cid: contour identifier or list of identifiers,

    :param str fname: output filename.

    :returns: None

    All contour segments will be saved to a zone called "Contour".
    Additional zones will be created for all segments with same
    boundary type.
    """
    icheck(0, UListOr1(ACont2D()))
    icheck(1, String())

    cb = flow.interface.ask_for_callback()
    cont = _cont2_from_id(cid)
    bt = flow.receiver.get_zone_types()
    tecplot_export.cont2(fname, cont, bt, cb)


@hmscriptfun
def export3d_surface_hmc(sid, fname, fmt="ascii"):
    """Exports 3d surface to hybmesh native format.

      :param sid: single or list of surface identifiers

      :param str fname: output filename

      :param str fmt: output data format:

         * ``'ascii'`` - all fields will be saved as text fields,
         * ``'bin'`` - all fields will be saved in binary section,
         * ``'fbin'`` - only floating point fields will be saved
           in binary section.

      See :ref:`surface3d-file` for format description.
    """
    if fmt == "binary":
        fmt = "bin"
    elif fmt == "fbinary":
        fmt = "fbin"
    icheck(0, UListOr1(ASurf3D()))
    icheck(1, String())
    icheck(2, OneOf('ascii', 'bin', 'fbin'))

    names = sid if isinstance(sid, list) else [sid]
    surfs = [_surf3_from_id(n) for n in names]
    cb = flow.interface.ask_for_callback()
    native_export.surf3_tofile(fname, surfs, names, fmt, cb)


@hmscriptfun
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
    if fmt == "binary":
        fmt = "bin"
    elif fmt == "fbinary":
        fmt = "fbin"
    icheck(0, String())
    icheck(1, OneOf('ascii', 'bin', 'fbin'))

    cb = flow.interface.ask_for_callback()
    native_export.export_all_tofile(fname, flow.receiver, fmt, cb)


@hmscriptfun
def save_project(fname, fmt="ascii"):
    """Saves current command flow and data to HybMesh project file.

       :param str fname: file name

       :param str fmt: output data format:

         * ``'ascii'`` - all fields will be saved as text fields,
         * ``'bin'`` - all fields will be saved in binary section,
         * ``'fbin'`` - only floating point fields will be saved
           in binary section.

       See :ref:`hmp-file` for description.
    """
    if fmt == "binary":
        fmt = "bin"
    elif fmt == "fbinary":
        fmt = "fbin"
    icheck(0, String())
    icheck(1, OneOf('ascii', 'bin', 'fbin'))

    flow_export.flow_and_framework_tofile(fname, flow, fmt)
