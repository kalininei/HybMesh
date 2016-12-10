import hybmeshpack.basic.proc as bp
import btypes


class Framework(object):
    """ Presents data collection for single work flow.
        Implements com.command.CommandReceiver interface
    """

    def __init__(self):
        #grids list
        self.grids2 = bp.NamedList()
        #user contours list
        self.contours2 = bp.NamedList()
        #grids3d list
        self.grids3 = bp.NamedList()
        #surfaces3 list
        self.surfaces3 = bp.NamedList()
        #boundary types
        self.boundary_types = btypes.BndTypesList()

    #boundary types manager
    def get_bnd_types(self):
        '-> btypes.BndTypesList'
        return self.boundary_types

    #grids manager
    def add_grid(self, name, grid):
        "-> (index, name, grid). add a grid to a grid list"
        self.grids2[name] = grid
        return self.grids2.get(val=grid)

    def remove_grid(self, name):
        'removes a grid with the certain name from a grid list'
        self.grids2.pop(name)

    def get_grid(self, ind=None, name=None, grid=None):
        "-> (index, name, grid). Get a grid by name or reference"
        return self.grids2.get(ind, name, grid)

    def get_grid_names(self):
        '-> [list of grid names]'
        return self.grids2.keys()

    #grids3d manager
    def add_grid3(self, name, grid):
        "-> (index, name, grid). add a grid to a grid list"
        self.grids3[name] = grid
        return self.grids3.get(val=grid)

    def remove_grid3(self, name):
        'removes a grid with the certain name from a grid list'
        self.grids3.pop(name)

    def get_grid3(self, ind=None, name=None, grid=None):
        "-> (index, name, grid). Get a grid by name or reference"
        return self.grids3.get(ind, name, grid)

    def get_grid3_names(self):
        '-> [list of grid names]'
        return self.grids3.keys()

    #contours manager
    def add_ucontour(self, name, cont):
        "-> (index, name, contour). add a user contour"
        self.contours2[name] = cont
        return self.contours2.get(val=cont)

    def remove_ucontour(self, name):
        "removes contour from framework"
        del self.contours2[name]

    def get_ucontour(self, ind=None, name=None, cont=None):
        "-> (index, name, grid). Get a contour by name or reference"
        return self.contours2.get(ind, name, cont)

    def get_ucontour_names(self):
        '-> [list of user contour names]'
        return self.contours2.keys()

    #surfaces manager
    def add_usurface(self, name, srf):
        "-> (index, name, contour). add a user surface"
        self.surfaces3[name] = srf
        return self.surfaces3.get(val=srf)

    def remove_usurface(self, name):
        "removes surface from framework"
        del self.surfaces3[name]

    def get_usurface(self, ind=None, name=None, srf=None):
        "-> (index, name, srf). Get a surface by name or reference"
        return self.surfaces3.get(ind, name, srf)

    def get_usurface_names(self):
        '-> [list of user surface names]'
        return self.surfaces3.keys()

    #all 2d names
    def get_all_names2(self):
        '-> [list of all grid and user contours]'
        return self.get_ucontour_names() + self.get_grid_names()

    #all 3d names
    def get_all_names3(self):
        '-> [list of all grids, contours, surfaces in 3d]'
        return self.get_usurface_names() +\
            self.get_grid3_names()

    #grid + user contour
    def get_any_contour(self, name):
        '-> Grid contour or user contour by its name'
        if name in self.contours2:
            return self.contours2[name]
        elif name in self.grids2:
            return self.grids2[name].cont
        else:
            return None

    #grid3 + user surface
    def get_any_surface(self, name):
        '-> Grid surface or user contour by its name'
        if name in self.surfaces3:
            return self.surfaces3[name]
        elif name in self.grids3:
            return self.grids3[name].surface()
        else:
            return None

    def get_any(self, name):
        '-> any 2d object by name'
        if name in self.contours2:
            return self.contours2[name]
        elif name in self.grids2:
            return self.grids2[name]
        else:
            return None

    #overriden from CommandReceiver
    def to_zero_state(self):
        'deletes all grids and contours'
        self.grids2.clear()
        self.contours2.clear()
        self.grids3.clear()
        self.surfaces3.clear()
        self.boundary_types.clear()

    def deep_copy(self):
        ret = Framework()
        #boundarys
        ret.boundary_types.add_data(self.boundary_types)
        #grids
        for name, g in self.grids2.items():
            gnew = g.deepcopy()
            ret.grids2[name] = gnew
        #contours
        for name, g in self.contours2.items():
            gnew = g.deepcopy()
            ret.contours2[name] = gnew
        #grids3
        for name, g in self.grids3.items():
            gnew = g.deepcopy()
            ret.grids3[name] = gnew
        #surfaces3
        for name, g in self.surfaces3.items():
            gnew = g.deepcopy()
            ret.surfaces3[name] = gnew
        return ret
