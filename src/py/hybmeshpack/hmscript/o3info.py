"3D objects information"
from hybmeshpack.hmscript import flow, hmscriptfun
from datachecks import (icheck, Grid3D, ASurf3D, Surf3D, OneOf)


@hmscriptfun
def info_grid3d(gid):
    """ Get 3d grid structure information

    :param gid: 3d grid identifier

    :returns: dictionary which represents
       total number of nodes, edges, faces, cells::

         {'Nnodes': int,
          'Nedges': int,
          'Nfaces': int,
          'Ncells': int
          }

    """
    icheck(0, Grid3D())
    g = flow.receiver.get_grid3(gid)
    ret = {}
    ret['Nnodes'] = g.n_vertices()
    ret['Nedges'] = g.n_edges()
    ret['Nfaces'] = g.n_faces()
    ret['Ncells'] = g.n_cells()
    return ret


@hmscriptfun
def info_surface(sid):
    """Get surface structure information

    :param sid: surface or 3d grid identifier

    :returns:
       dictionary representing total number of nodes, edges, faces,
       and number of faces of each boundary type::

         {'Nnodes': int,
          'Nedges': int,
          'Nfaces': int,
          'btypes': {btype(int): int}  # boundary type: number of faces
         }
    """
    icheck(0, ASurf3D())

    s = flow.receiver.get_any_surface(sid)
    ret = {}
    ret['Nnodes'] = s.n_vertices()
    ret['Nedges'] = s.n_edges()
    ret['Nfaces'] = s.n_faces()
    ret['btypes'] = {}
    for b in s.raw_data('bt'):
        if b not in ret['btypes']:
            ret['btypes'][b] = 0
        ret['btypes'][b] += 1
    return ret


@hmscriptfun
def domain_volume(sid):
    """Calculates area of closed domain bounded by the given surface

    :param sid: grid3d or surface identifier

    :returns: positive float or zero for not closed surfaces
    """
    icheck(0, ASurf3D())
    return flow.receiver.get_any_surface(sid).volume()


@hmscriptfun
def tab_surf3(obj, what):
    """ Returns plain table for the given surface.

    :param str obj: surface identifier

    :param str what: table name

    :returns: plain ctypes array representing requested table.

    Possible **what** values are:

    * ``'vert'`` - vertex coordinates table,
    * ``'face_center'`` - faces center point coordinates table,
    * ``'edge_vert'`` - edge-vertex connectivity: indices of first and
      last vertices for each edge,
    * ``'face_dim'`` - number of vertices in each face,
    * ``'face_edge'`` - face-edge connectivity: ordered list of edge indicies
      for each face
    * ``'face_vert'`` - face-vertex connectivity: ordered list of vertex
      indices for each face
    * ``'bt'`` - boundary features of each edge.

    In case of surfaces with variable face dimensions
    ``'face_edge'`` and ``'face_vert'`` tables require
    additional ``'face_dim'`` table to subdive
    returned plain array by certain faces.

    """
    icheck(0, Surf3D())
    icheck(1, OneOf('vert',
                    'edge_vert',
                    'face_dim', 'face_edge', 'face_vert',
                    'bt',
                    'face_center'
                    ))
    return flow.receiver.get_surface3(obj).raw_data(what)


@hmscriptfun
def tab_grid3(obj, what):
    """ Returns plain table for the given grid.

    :param str obj: grid identifier

    :param str what: table name

    :returns: plain ctypes array representing requested table.

    Possible **what** values are:

    * ``'vert'`` - vertex coordinates table,
    * ``'edge_vert'`` - edge-vertex connectivity: indices of first and
      last vertices for each edge,
    * ``'face_dim'`` - number of vertices in each face,
    * ``'face_edge'`` - face-edge connectivity: ordered list of edge indicies
      for each face
    * ``'face_vert'`` - face-vertex connectivity: ordered list of vertex
      indicies for each face
    * ``'face_cell'`` - face-cell connectivity: indices of left and right
      cell for each face. If this is a boundary face ``-1`` is used to mark
      boundary side,
    * ``'cell_fdim'`` - number of faces in each cell,
    * ``'cell_vdim'`` - number of vertices in each cell,
    * ``'cell_face'`` - cell-edge connectivity: edge indices for each cell.
    * ``'cell_vert'`` - cell-vertex connectivity: vertex indices for each cell.
      Vertices of simple shaped cells are ordered according to vtk
      file format. See figure in :func:`export3d_grid_hmg`.
    * ``'bnd'`` - list of boundary faces indices,
    * ``'bt'`` - boundary types for all faces including internal ones,
    * ``'bnd_bt'`` - (boundary face, boundary feature) pairs

    Use ``'face_dim'``, ``'cell_fdim'``, ``'cell_vdim'`` to subdivide
    plain ``'face_vert'``, ``'face_edge'``, ``'cell_face'``, ``'cell_vert'``
    arrays by certain faces/cells.
    """
    icheck(0, Grid3D())
    icheck(1, OneOf('vert',
                    'edge_vert',
                    'face_dim', 'face_edge', 'face_vert', 'face_cell',
                    'cell_fdim', 'cell_vdim', 'cell_face', 'cell_vert',
                    'bt', 'bnd', 'bnd_bt',
                    ))
    return flow.receiver.get_grid3(obj).raw_data(what)
