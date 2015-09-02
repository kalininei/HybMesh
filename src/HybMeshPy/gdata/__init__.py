import basic.proc as bp
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
        #boundary types
        self.boundary_types = btypes.BndTypesList()

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

    #boundary types manager
    def get_bnd_types(self):
        '-> btypes.BndTypesList'
        return self.boundary_types

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

    #grid + contours
    def get_any_contour(self, name):
        '-> Grid or user contour by its name'
        if name in self.contours2:
            return self.contours2[name]
        elif name in self.grids2:
            return self.grids2[name].cont
        else:
            return None

    def get_all_names(self):
        '-> [list of all grid and user contours]'
        return self.get_contour_names() + self.get_grid_names()

    #overriden from CommandReceiver
    def to_zero_state(self):
        'deletes all grids and contours'
        self.grids2.clear()
        self.contours2.clear()
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
        return ret
