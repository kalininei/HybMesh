import command
import objcom
import gridcom
import contcom
# import importcom


class _Factory(object):
    ' Command class Factory '
    command_list = [
            command.Start,
            gridcom.AddUnfRectGrid,
            gridcom.AddUnfCircGrid,
            gridcom.AddUnfRingGrid,
            gridcom.UniteGrids,
            gridcom.ExcludeContours,
            # gridcom.BuildBoundaryGrid,
            objcom.RenameGeom,
            objcom.RemoveGeom,
            objcom.MoveGeom,
            objcom.RotateGeom,
            objcom.ScaleGeom,
            objcom.CopyGeom,
            contcom.AddRectCont,
            contcom.UniteContours,
            contcom.EditBoundaryType,
            contcom.GridBndToContour,
            contcom.SimplifyContours,
            contcom.SetBTypeToContour,
            # importcom.ImportContourNative,
            # importcom.ImportContourASCII,
            # importcom.ImportGridNative,
            # importcom.ImportGridSimple34,
                   ]

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
        self._coms[comcls.method_code()] = comcls

factory = _Factory()
