import command
import objcom
import gridcom
import grid3dcom
import contcom
import imcom


def _get_classes(module):
    import inspect
    ret = []
    for name, obj in inspect.getmembers(module):
        if inspect.isclass(obj) and obj.__name__[0] != '_':
            ret.append(obj)
    return ret


class _Factory(object):
    ' Command class Factory '
    command_list = [command.Start] + _get_classes(objcom) +\
        _get_classes(gridcom) + _get_classes(grid3dcom) +\
        _get_classes(contcom) + _get_classes(imcom)

    def __init__(self):
        """ Initializes the commands building factory.
            coms - is the list of classes which
            implement Command interface
        """
        self._coms = {}
        for c in self.command_list:
            self.__register(c)

    def string_from_cls(self, cls):
        ind = self._coms.values().index(cls)
        return self._coms.keys()[ind]

    def cls_from_string(self, s):
        return self._coms[s]

    def create_from_string(self, title, strcom):
        cls = self.cls_from_string(title)
        return cls.fromstring(strcom)

    def create_from_args(self, clscode, *args):
        cls = self.cls_from_string(clscode)
        return cls(*args)

    def __register(self, comcls):
        try:
            self._coms[comcls.method_code()] = comcls
        except:
            pass

factory = _Factory()
