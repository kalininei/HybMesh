import objcom
import command


class Grid3BndToSurface(objcom.AbstractAddRemove):
    'Copy grid boundary to user surface structure'

    def __init__(self, arg):
        if "surf_name" not in arg:
            arg["surf_name"] = "Surface1"
        if "separate" not in arg:
            arg['separate'] = False
        super(Grid3BndToSurface, self).__init__(arg)

    def doc(self):
        return "Copy boundary from %s" % self.options['grid_name']

    @classmethod
    def _arguments_types(cls):
        return {'grid_name': command.BasicOption(str),
                'surf_name': command.BasicOption(str),
                'separate': command.BoolOption(),
                }

    def _addrem_objects(self):
        try:
            g = self.grid3_by_name(self.options['grid_name'])
            s = [g.surface().deepcopy()]
            ret = []
            if self.options['separate']:
                s = s[0].shallow_separate()
            for it in s:
                ret.append((self.options['surf_name'], it))
            return [], [], [], [], [], [], ret, []
        except Exception as e:
            raise command.ExecutionError(str(e), self, e)
