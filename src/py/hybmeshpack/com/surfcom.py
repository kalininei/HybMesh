import addremove
import comopt as co
from hybmeshpack.hmcore import s3 as s3core
from hybmeshpack.gdata.srf3 import Surface3


class Grid3BndToSurface(addremove.AbstractAddRemove):
    'Copy grid boundary to user surface structure'

    def __init__(self, arg):
        super(Grid3BndToSurface, self).__init__(arg)

    def doc(self):
        return "Copy boundary from %s" % self.options['grid_name']

    @classmethod
    def _arguments_types(cls):
        return {'surf_name': co.BasicOption(str, None),
                'grid_name': co.BasicOption(str),
                'separate': co.BoolOption(False),
                }

    def _addrem_surface3(self):
        gn = self.get_option('grid_name')
        sn = self.get_option('surf_name')
        if sn is None:
            sn = gn + '_surface'
        g = self.grid3_by_name(gn)
        s = [g.surface().surface3()]
        if self.get_option('separate'):
            s = s3core.quick_separate(s[0].cdata)
            s = [Surface3(it) for it in s]

        return [(it, sn) for it in s], []
