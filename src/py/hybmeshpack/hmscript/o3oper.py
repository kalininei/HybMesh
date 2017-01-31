" 3D objects operations"
from hybmeshpack import com
from hybmeshpack.hmscript import flow, ExecError


def merge_grids3(g1, g2):
    """ Merges 3d grids into single one.

        :param g1:

        :param g2: 3d source grids identifiers.

        :returns: new grid identifier.

        :raises: ExecError

        Merge procedure will process only strictly
        coincident boundary primitives.
    """
    c = com.grid3dcom.Merge({"src1": g1, "src2": g2})
    try:
        flow.exec_command(c)
        return c.added_grids3()[0]
    except:
        raise ExecError('merge_grids3')
