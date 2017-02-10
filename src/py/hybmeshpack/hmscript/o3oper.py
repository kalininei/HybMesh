" 3D objects operations"
from hybmeshpack import com
from hybmeshpack.hmscript import flow, hmscriptfun
from datachecks import icheck, Grid3D


@hmscriptfun
def merge_grids3(g1, g2):
    """ Merges 3d grids into single one.

        :param g1:

        :param g2: 3d source grids identifiers.

        :returns: new grid identifier.

        Merge procedure will process only strictly
        coincident boundary primitives.
    """
    icheck(0, Grid3D())
    icheck(1, Grid3D())

    c = com.grid3dcom.Merge({"src1": g1, "src2": g2})
    flow.exec_command(c)
    return c.added_grids3()[0]
