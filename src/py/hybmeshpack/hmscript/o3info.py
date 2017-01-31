"3D objects information"
from hybmeshpack.hmscript import flow


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
    g = flow.receiver.get_grid3(gid)
    ret = {}
    ret['Nnodes'] = g.n_points()
    ret['Nedges'] = g.n_edges()
    ret['Nfaces'] = g.n_faces()
    ret['Ncells'] = g.n_cells()
    return ret


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

    s = flow.receiver.get_any_surface(sid)
    ret = {}
    ret['Nnodes'] = s.n_points()
    ret['Nedges'] = s.n_edges()
    ret['Nfaces'] = s.n_faces()
    ret['btypes'] = {}
    for b in s.raw_data('btypes'):
        if b not in ret['btypes']:
            ret['btypes'][b] = 0
        ret['btypes'][b] += 1
    return ret


def domain_volume(sid):
    """Calculates area of closed domain bounded by the **s** surface

    :param sid: grid3d or surface identifier

    :returns: positive float or zero for not closed surfaces
    """
    return flow.receiver().get_any_surface(sid).volume()
