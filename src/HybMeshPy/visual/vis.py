import globvars
import contvis
import gridvis
import btypes
import geomdlgs


class FrameworkVis(object):
    'Interface for visual representation of framework object'
    def __init__(self):
        self.fw = None

    def set_framework(self, fw):
        'sets a reference to the framework object @fw '
        self.fw = fw

    def update(self):
        ' Updates view of referenced framework '
        raise NotImplementedError

    #requests
    def ask_checked_grid_names(self):
        'returns names of checked grids '
        raise NotImplementedError

    def ask_checked_contours_names(self):
        'returns names of checked contours'
        raise NotImplementedError

    def ask_checked_gbnd_names(self):
        'returns names of checkes grid boundaries'
        raise NotImplementedError

    def ask_new_contours_bnd(self, cnames, btypes):
        """ -> [{boundry-index: [list of contour edges]}, {}, ...]
            -> None
            Request for new boundary types for contours.
        """
        raise NotImplementedError

    #checks and visibility
    def check_grid(self, name, yesno):
        """ Checks (yesno=True) or unchecks (False) the grid by its name.
            Check/Uncheck all grids if name is None
        """
        raise NotImplementedError

    def check_cont(self, name, yesno):
        """ Checks (yesno=True) or unchecks (False) the contour by its name.
            Check/Uncheck all contours if name is None
        """
        raise NotImplementedError

    def show_grid(self, name, yesno):
        """ Shows (yesno=True) or Hide (False) the grid by its name.
            Show/Hide all grids if name is None
        """
        raise NotImplementedError

    def show_cont(self, name, yesno):
        """ Shows (yesno=True) or Hide (False) the contour by its name.
            Show/Hide all contours if name is None
        """
        raise NotImplementedError


class QtFrameworkVis(FrameworkVis):
    """ output into three Qt Frames:
         -- visual
         -- grid tree
         -- contours tree
    """
    def __init__(self):
        super(QtFrameworkVis, self).__init__()
        self.cont_view = globvars.mainWindow.cont_view
        self.grid_view = globvars.mainWindow.grid_view
        self.bnd_view = globvars.mainWindow.bnd_view
        self.grid_model = None
        self.cont_model = None
        self.bnd_model = None
        self.vtk_draw = globvars.mainWindow.vtkWidget

    #======================== OVERRIDEN
    def set_framework(self, fw):
        super(QtFrameworkVis, self).set_framework(fw)
        self.cont_model = contvis.ContListModel(self, self.cont_view)
        self.grid_model = gridvis.GridListModel(self, self.grid_view)
        self.bnd_model = btypes.BoundaryTypesModel(fw.boundary_types,
                self.bnd_view)

    def ask_checked_grid_names(self):
        ' returns names of checked grids '
        try:
            return self.grid_model.checked_names()
        except:
            return []

    def ask_checked_contours_names(self):
        'returns names of checked contours'
        try:
            return self.cont_model.checked_names_udc()
        except:
            return []

    def ask_checked_gbnd_names(self):
        'returns names of checkes grid boundaries'
        try:
            return self.cont_model.checked_names_gbc()
        except:
            return []

    def ask_new_contours_bnd(self, cnames, btypes):
        """ ->[{boundry-index: [list of contour edges]}, {}, ...]
            -> None
            Request for new boundary types for contours.
        """
        acts = [self.cont_model.get_actor(n) for n in cnames]
        polys = [a.GetMapper().GetInput() for a in acts]
        conts = [self.fw.get_any_contour(n) for n in cnames]
        arg = zip(cnames, conts, polys)
        dialog = geomdlgs.BoundaryTypesManagement(arg, btypes)
        if dialog.exec_():
            return dialog.ret_value()
        else:
            return None

    def check_grid(self, name, yesno):
        g = self.fw.get_grid(name=name)[2] if name is not None else None
        self.grid_model.check(g, yesno)

    def check_cont(self, name, yesno):
        g = self.fw.get_any_contour(name) if name is not None else None
        self.cont_model.check(g, yesno)

    def show_grid(self, name, yesno):
        g = self.fw.get_grid(name=name)[2] if name is not None else None
        self.grid_model.show(g, yesno)

    def show_cont(self, name, yesno):
        g = self.fw.get_any_contour(name) if name is not None else None
        self.cont_model.show(g, yesno)

    def update(self):
        """ updates view in of referenced framework in the main window '
            after new command has been executed
        """
        #cold start invocation before window.show() call
        if not globvars.mainWindow.isVisible():
            return

        #============== Grids
        if self.grid_view.model() is not self.grid_model:
            self.grid_view.setModel(self.grid_model)
        self.grid_model.reset()

        #============== CONTOURS
        if self.cont_view.model() is not self.cont_model:
            self.cont_view.setModel(self.cont_model)
        self.cont_model.reset()
        self.cont_view.expandAll()

        #============== Boundary types
        if self.bnd_view.model() is not self.bnd_model:
            self.bnd_view.setModel(self.bnd_model)
        self.bnd_model.reset()

        #============== Draw
        #set actors list
        alist = self.grid_model.active_actors() +\
                    self.cont_model.active_actors()
        self.vtk_draw.set_actor_list(alist)
        self.vtk_draw.ren.ResetCamera()
        self.vtk_draw.Render()
