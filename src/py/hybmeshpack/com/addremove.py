import command


class _AddRemoveObjects(object):
    """ subcommand to remove objects.
        Could be used only as the internal command for other user commands

        Order of execution:
            remove grids2, contours2, grids3, surface3
            add grids2, contours2, grids3, surface3
    """
    def __init__(self,
                 addg2=[], remg2=[],
                 addc2=[], remc2=[],
                 addg3=[], remg3=[],
                 adds3=[], rems3=[]):
        """ addg2 -- [ (Grid2, name), ... ]
            remg2 -- [ string name, ... ]
            addc2 -- [ (Contour2, name), ...]
            remc2 -- [ string name, ...]
            addg3    -- [ (Grid3, name), ... ]
            remg3    -- [ string name, ...]
            adds3    -- [ (Surface srf, name), ... ]
            rems3    -- [ string name, ...]
        """
        self.clear_backup()
        self.addg2, self.addc2 = addg2, addc2
        self.remg2, self.remc2 = remg2, remc2
        self.addg3, self.remg3 = addg3, remg3
        self.adds3, self.rems3 = adds3, rems3

    def do(self, receiver):
        self.receiver = receiver
        if len(self.addg2) + len(self.remg2) +\
                len(self.addc2) + len(self.remc2) +\
                len(self.addg3) + len(self.remg3) +\
                len(self.adds3) + len(self.rems3) == 0:
            return False
        self.clear_backup()
        self.backup = receiver.backup_copy()

        # remove grids
        for n in self.remg2:
            receiver.remove_grid2(name=n)
        # remove contours
        for n in self.remc2:
            receiver.remove_contour2(name=n)
        # remove g3
        for n in self.remg3:
            receiver.remove_grid3(name=n)
        # remove s3
        for n in self.rems3:
            receiver.remove_surface3(name=n)
        # add grids
        for v in self.addg2:
            self.new_g2.append(receiver.append_grid2(*v))
        # add contours
        for v in self.addc2:
            self.new_c2.append(receiver.append_contour2(*v))
        # add grid3
        for v in self.addg3:
            self.new_g3.append(receiver.append_grid3(*v))
        # add surface3
        for v in self.adds3:
            self.new_s3.append(receiver.append_surface3(*v))
        return True

    def undo(self):
        self.receiver.restore_from_backup(self.backup)
        self.clear_backup()

    def clear_backup(self):
        self.new_g2 = []
        self.new_g3 = []
        self.new_c2 = []
        self.new_s3 = []
        self.backup = None


class AbstractAddRemove(command.Command):
    """ Abstract base for commands which end up
        with adding and/or removing objects """
    def __init__(self, argsdict):
        super(AbstractAddRemove, self).__init__(argsdict)

    def _exec(self):
        self.__addrem = _AddRemoveObjects(*self._addrem_objects())
        return self.__addrem.do(self.receiver)

    def _clear(self):
        self.__addrem = None

    def _undo(self):
        self.__addrem.undo()

    def _redo(self):
        self.__addrem.do(self.receiver)

    def added_grids2(self):
        '->list of new 2d grids names'
        return self.__addrem.new_g2

    def added_contours2(self):
        '->list of new 2d contours names'
        return self.__addrem.new_c2

    def added_grids3(self):
        '->list of new 3d grids names'
        return self.__addrem.new_g3

    def added_surfaces3(self):
        '->list of new 3d surfaces names'
        return self.__addrem.new_s3

    def _addrem_objects(self):
        """ -> addg2, remg2, addc2, remc2, addg3, remg3, adds3, rems3
        addg2 -- [ (Grid2, name), ... ]
        remg2 -- [ string name, ... ]
        addc2 -- [ (Contour2, name), ...]
        remc2 -- [ string name, ...]
        addg3    -- [ (Grid3, name), ... ]
        remg3    -- [ string name, ...]
        adds3    -- [ (Surface srf, name), ... ]
        rems3    -- [ string name, ...]
        """
        return self._addrem_grid2() + self._addrem_contour2() +\
            self._addrem_grid3() + self._addrem_surface3()

    # function for overriding
    def _addrem_grid2(self):
        """ -> addg2, remg2
        addg2 -- [ (Grid2, name), ... ]
        remg2 -- [ string name, ... ]
        """
        return [], []

    def _addrem_contour2(self):
        """ -> addc2, remc2
        addc2 -- [ (Contour2, name), ... ]
        remc2 -- [ string name, ... ]
        """
        return [], []

    def _addrem_grid3(self):
        """ -> addg3, remg3
        addg3 -- [ (Grid3, name), ... ]
        remg3 -- [ string name, ... ]
        """
        return [], []

    def _addrem_surface3(self):
        """ -> adds3, rems3
        adds3 -- [ (Surface3, name), ... ]
        rems3 -- [ string name, ... ]
        """
        return [], []
