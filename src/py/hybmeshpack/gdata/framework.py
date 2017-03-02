import itertools
from basic import BndType
from hybmeshpack.basic.namedlist import NamedList, NamedListPool
from hybmeshpack.gdata import grid2, cont2, srf3, grid3


class Framework(object):
    """ Data collection for single work flow  """

    def __init__(self):
        self._obj_pool = NamedListPool()
        self._znames_pool = NamedListPool()
        # grids list
        self.grids2 = NamedList("Grid2D_1", self._obj_pool)
        # user contours list
        self.contours2 = NamedList("Contour2D_1", self._obj_pool)
        # grids3d list
        self.grids3 = NamedList("Grid3D_1", self._obj_pool)
        # surfaces3 list
        self.surfaces3 = NamedList("Surface3D_1", self._obj_pool)
        # zone types
        self.zone_types = NamedList("Zone1", self._znames_pool)
        self.zone_types.append(BndType(0), "default-boundary")

    # boundary types manager
    def get_zone_types(self):
        '-> {zone type: zone name}'
        ret = {}
        for n in self.zone_types.all_names():
            ind = self.zone_types.obj_by_name(n).index
            ret[ind] = n
        return ret

    def remove_zone_type(self, remi):
        if remi == 0:
            return
        for ob in self.zone_types.all_obj():
            if ob.index == remi:
                self.zone_types.remove(obj=ob)
                return

    def set_zone_type(self, index, name=None):
        if name is None:
            name = 'boundary' + str(index)
        for ob in self.zone_types.all_obj():
            if ob.index == index:
                self.zone_types.change(oldobj=ob, newname=name)
                return
        self.zone_types.append(BndType(index), name)

    # 2d grids manager
    def append_grid2(self, grid, name=None):
        "-> name. add a grid to a grid list"
        return self.grids2.append(grid, name)

    def remove_grid2(self, grid=None, name=None):
        'removes a grid with the certain name from a grid list'
        self.grids2.remove(grid, name)

    def get_grid2(self, name):
        "-> grid. Get a grid by name"
        return self.grids2.obj_by_name(name)

    def get_grid2_names(self):
        '-> [list of grid names]'
        return self.grids2.all_names()

    # 3d grids manager
    def append_grid3(self, grid, name=None):
        "-> name. add a grid to a grid list"
        return self.grids3.append(grid, name)

    def remove_grid3(self, grid=None, name=None):
        'removes a grid with the certain name from a grid list'
        self.grids3.remove(grid, name)

    def get_grid3(self, name):
        "-> grid. Get a grid by name"
        return self.grids3.obj_by_name(name)

    def get_grid3_names(self):
        '-> [list of grid names]'
        return self.grids3.all_names()

    # 2d contours manager
    def append_contour2(self, contour, name=None):
        "-> name. add a contour to the contour list"
        return self.contours2.append(contour, name)

    def remove_contour2(self, contour=None, name=None):
        'removes a contour with the certain name from the contour list'
        self.contours2.remove(contour, name)

    def get_contour2(self, name):
        "-> contour. Get a contour by name"
        return self.contours2.obj_by_name(name)

    def get_contour2_names(self):
        '-> [list of contour names]'
        return self.contours2.all_names()

    # 3d surfaces manager
    def append_surface3(self, surface, name=None):
        "-> name. add a surface to the surface list"
        return self.surfaces3.append(surface, name)

    def remove_surface3(self, surface=None, name=None):
        'removes a surface with the certain name from the surface list'
        self.surfaces3.remove(surface, name)

    def get_surface3(self, name):
        "-> surface. Get a surface by name or reference"
        return self.surfaces3.obj_by_name(name)

    def get_surface3_names(self):
        '-> [list of surface names]'
        return self.surfaces3.all_names()

    # 2d objects manager
    def get_names2(self):
        '-> [list of all grid and contours names]'
        return self.get_grid2_names() + self.get_contour2_names()

    def get_object2(self, name):
        '-> GeomObject2 by its name'
        try:
            return self.grids2.obj_by_name(name)
        except KeyError:
            return self.contours2.obj_by_name(name)

    def get_any_contour(self, name):
        '-> AbstractContour2. Returns contour or grid contour object'
        try:
            return self.grids2.obj_by_name(name).contour()
        except KeyError:
            return self.contours2.obj_by_name(name)

    # 3d objects manager
    def get_names3(self):
        '-> [list of all grid and surfaces names]'
        return self.get_grid3_names() + self.get_surface3_names()

    def get_object3(self, name):
        '-> GeomObject3 by its name'
        try:
            return self.grids3.obj_by_name(name)
        except KeyError:
            return self.surfaces3.obj_by_name(name)

    def get_any_surface(self, name):
        '-> AbstractSurface. Returns surface or grid surface object'
        try:
            return self.grids3.obj_by_name(name).surface()
        except KeyError:
            return self.surfaces3.obj_by_name(name)

    # all object manager
    def get_names(self):
        '-> [list of all names]'
        return self.get_grid2_names() + self.get_grid3_names() +\
            self.get_contour2_names() + self.get_surface3_names()

    def get_object(self, name):
        '-> GeomObject by name'
        try:
            return self.get_object2(name)
        except KeyError:
            return self.get_object3(name)

    def whatis(self, iden):
        if isinstance(iden, str):
            iden = self.get_object(iden)
        if isinstance(iden, cont2.Contour2):
            return 'c2'
        elif isinstance(iden, grid2.Grid2):
            return 'g2'
        elif isinstance(iden, srf3.Surface3):
            return 's3'
        elif isinstance(iden, grid3.Grid3):
            return 'g3'
        elif isinstance(iden, grid3.GridSurface):
            return 'g3s'
        elif isinstance(iden, grid2.GridContour):
            return 'g2c'
        else:
            raise KeyError

    # state manipulations
    def to_zero_state(self):
        'deletes everything'
        self.grids2.clear()
        self.contours2.clear()
        self.grids3.clear()
        self.surfaces3.clear()
        self.zone_types.clear()
        self._znames_pool.clear()
        self._obj_pool.clear()
        self.set_zone_type(0, "default-boundary")

    # copy
    def shallow_fillfrom(self, fwork):
        """ get data from another framework.
        Uses shallow copies of objects and deep copies for all names
        and dictionaries.
        """
        self._obj_pool.fillfrom(fwork._obj_pool)
        self._znames_pool.fillfrom(fwork._znames_pool)
        self.grids2.shallow_fillfrom(fwork.grids2)
        self.grids3.shallow_fillfrom(fwork.grids3)
        self.contours2.shallow_fillfrom(fwork.contours2)
        self.surfaces3.shallow_fillfrom(fwork.surfaces3)
        self.zone_types.shallow_fillfrom(fwork.zone_types)

    def backup_copy(self):
        """ -> Framework which could be used for state restoration.
        It shares all objects whith this one but has its own deep copies
        of NamedLists and NamedListPools.
        You could safely add and remove new objects in one not disturbing
        the other. However changes of objects affect both frameworks.
        """
        ret = Framework()
        ret.shallow_fillfrom(self)
        return ret

    def restore_from_backup(self, backup):
        """ restores state from backup copy """
        self.shallow_fillfrom(backup)

    def deepcopy(self):
        """ -> Framework. Makes deep copy of all objects """
        ret = Framework()
        ret.shallow_fillfrom(self)
        ch = [ret.grids2._data, ret.grids3._data, ret.contours2._data,
              ret.surfaces3._data, ret.zone_types._data]
        for e in itertools.chain(*ch):
            e.obj = e.obj.deepcopy()
        return ret
