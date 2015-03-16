import sip

from PyQt4 import QtGui
import vtk

import qtbp
import gridcom
import objcom
import globvars
import dlgs
import contvis
import btypes
import geomdlgs

#vtk options
def_color = (0.7, 0.7, 0.7)
def_linewidth = 1
chk_color = (1, 0, 0)
chk_linewidth = 2


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

    #user requests
    def ask_checked_grid_names(self):
        ' returns names of checked grids '
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


class GridBlockTreeItem(QtGui.QTreeWidgetItem):
    ' presents a gridblock item in a grid manager docker window '

    def __init__(self, name, view_opt):
        super(GridBlockTreeItem, self).__init__()
        self.view_opt = view_opt
        self.name = name
        #options
        self.setIcon(3, qtbp.get_icon("opts"))
        #delete button
        self.setIcon(4, qtbp.get_icon("del"))
        #opts
        self.update(name)

    def update(self, name=None):
        ' change view according to view_options '
        #text
        if name is not None:
            self.setText(0, name)
        #visibility
        self._set_view_icon()
        #check status
        if self.view_opt.get_checked():
            chs = 2
        else:
            chs = 0
        self.setCheckState(1, chs)

    def _set_view_icon(self):
        ' sets visibility icon '
        if self.view_opt.get_visibility():
            self.setIcon(2, qtbp.get_icon("eye-on"))
        else:
            self.setIcon(2, qtbp.get_icon("eye-off"))

    def clicked(self, col):
        """ item click handler
            col -- index of column
        """
        if col == 2:
            #visibility change
            self.view_opt.set_visibility(not self.view_opt.get_visibility())
            globvars.mainWindow.vtkWidget.GetRenderWindow().Render()
        elif col == 1:
            #check state change
            if self.checkState(col) == 2:
                chs = True
            else:
                chs = False
            self.view_opt.set_checked(chs)
        elif col == 3:
            #view options
            dlgs.GridViewOpt(None, globvars.mainWindow).exec_()
        elif col == 4:
            #delete grid
            #assuming that we are operating on the actual flow
            com = objcom.RemoveGrid([self.name])
            globvars.actual_flow().exec_command(com)

    def show_context(self, pnt):
        def trig(m, name, fun):
            "create menu item"
            act = QtGui.QAction(name, self.treeWidget())
            act.triggered.connect(fun)
            m.addAction(act)
        menu = QtGui.QMenu(self.treeWidget())
        trig(menu, "Rename", self._rename)
        menu.popup(self.treeWidget().viewport().mapToGlobal(pnt))

    def _rename(self):
        oldname = self.name
        newname, ok = QtGui.QInputDialog.getText(self.treeWidget(),
                "Rename grid", "Enter new grid name", text=oldname)
        oldname, newname = str(oldname), str(newname)
        if ok and oldname != newname:
            com = gridcom.RenameGrid(oldname, newname)
            globvars.actual_flow().exec_command(com)


class ViewGridData(object):
    """ A container for grid block representation
    Holds actors, Grid Manager QTreeWidgetItem

    All objects are unique for each grid_id and
    don't free memory with the deletion of the corresponding grid
    in order to support undo/redo strategy.
    Hence corresponding actors and other memory consuming attrs
    should be freed manualy with the deletion of the grid by
    invocation of mem_minimize() method
    """
    #id -> ViewGridData
    stored = {}

    def __new__(cls, g2):
        """ Creates new object if its id not in cls.stored.
            Otherwise returns object from stored[id]
        """
        if g2.get_id() in cls.stored:
            ret = cls.stored[g2.get_id()]
        else:
            ret = super(ViewGridData, cls).__new__(cls)
            #default options
            ret._visible = True
            ret._checked = False
            #data
            ret._grid_actor = None
            ret._manager_item = None
            #add to dictionary
            cls.stored[g2.get_id()] = ret

        #adds a subscriber if grid geometry was chanded
        g2.add_subscriber_change_geom(ret.grid_geom_event)
        return ret

    def mem_minimize(self):
        " frees all memory consuming data "
        self._grid_actor = None
        self._manager_item = None

    def grid_geom_event(self):
        """ if grid geometry was change -> delete an actor
            It will be created again in a get_actor method
        """
        self._grid_actor = None

    #---------------------- Options operations
    # ---- Visibility
    def get_visibility(self):
        return self._visible

    def set_visibility(self, isvis=None):
        if isvis is not None:
            self._visible = isvis
        if self._visible:
            self._grid_actor.VisibilityOn()
        else:
            self._grid_actor.VisibilityOff()
        self._manager_item.update()

    # ---- Check Status
    def get_checked(self):
        return self._checked

    def set_checked(self, chs=None):
        if chs is not None:
            self._checked = chs
        self._manager_item.update()
        if self._checked:
            self._grid_actor.GetProperty().SetColor(*chk_color)
            self._grid_actor.GetProperty().SetLineWidth(chk_linewidth)
        else:
            self._grid_actor.GetProperty().SetColor(*def_color)
            self._grid_actor.GetProperty().SetLineWidth(def_linewidth)
        globvars.mainWindow.vtkWidget.GetRenderWindow().Render()

    # --------------------- Data operations
    # ---- actor procedures
    def get_actor(self, g2):
        if self._grid_actor is None:
            self._reset_actor(g2)
        return self._grid_actor

    def _reset_actor(self, g2):
        " Create new actor for the grid "
        ps = vtk.vtkPoints()
        g = vtk.vtkUnstructuredGrid()
        g.Allocate(1, 1)  # number of cells?
        #points
        ps.SetNumberOfPoints(g2.n_points())
        for i, pnt in enumerate(g2.points):
            ps.InsertPoint(i, pnt.x, pnt.y, 0)
        g.SetPoints(ps)
        #cells
        for i, cl in enumerate(g2.cells_nodes_connect()):
            cell = vtk.vtkPolygon()
            cell.GetPointIds().Allocate(len(cl))
            #WTF? without reversing some edges are black
            map(cell.GetPointIds().InsertNextId, cl[::-1])
            g.InsertNextCell(cell.GetCellType(), cell.GetPointIds())

        #mapper and actor
        mp = vtk.vtkDataSetMapper()
        mp.SetInputData(g)
        self._grid_actor = vtk.vtkActor()
        self._grid_actor.GetProperty().SetRepresentationToWireframe()
        self._grid_actor.SetMapper(mp)

        #apply old visibility and checked options to new actor
        self.set_visibility()
        self.set_checked()

    # ----- Tree Item procedures
    def get_manager_item(self, name, grid):
        #since treeWidgetItem is freed from c++ memory
        #after removing from its tree we should also
        #check for existance of its c++ core using sip.isdeleted
        if self._manager_item is None or sip.isdeleted(self._manager_item):
            self._manager_item = GridBlockTreeItem(name, self)
        self._manager_item.update(name)
        return self._manager_item


class QtFrameworkVis(FrameworkVis):
    """ output into three Qt Frames:
         -- visual
         -- grid tree
         -- contours tree
    """
    def __init__(self):
        super(QtFrameworkVis, self).__init__()
        self.cont_view = globvars.mainWindow.cont_view
        self.bnd_view = globvars.mainWindow.bnd_view
        self.cont_model = None
        self.bnd_model = None
        self.vtk_draw = globvars.mainWindow.vtkWidget

    #======================== OVERRIDEN
    def set_framework(self, fw):
        super(QtFrameworkVis, self).set_framework(fw)
        self.cont_model = contvis.ContListModel(fw.contours2,
                fw.grids2, fw.boundary_types, self.cont_view)
        self.bnd_model = btypes.BoundaryTypesModel(fw.boundary_types,
                self.bnd_view)

    def ask_checked_grid_names(self):
        ' returns names of checked grids '
        ret = []
        for (name, grid) in self.fw.grids2.items():
            view_data = ViewGridData(grid)
            if view_data.get_checked():
                ret.append(name)
        return ret

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

    def update(self):
        """ updates view in of referenced framework in the main window '
            after new command has been executed
        """
        #cold start invocation before window.show() call
        if not globvars.mainWindow.isVisible():
            return

        #============== Grids
        #clear vtk actors for grids which are deleted
        ids = [a.get_id() for a in self.fw.grids2.values()]
        for k, v in ViewGridData.stored.items():
            if k not in ids:
                v.mem_minimize()
        #init list of active grid actors
        grid_actors = []
        #clear all items in grids manager window
        globvars.mainWindow.clear_grid_manager()
        for (name, grid) in self.fw.grids2.items():
            #create or restore view options for the grid
            view_data = ViewGridData(grid)
            #create the grid item for the grid manager window
            tree_item = view_data.get_manager_item(name, grid)
            #add the created item
            globvars.mainWindow.add_grid_manager_item(tree_item)
            #grid actors list
            grid_actors.append(view_data.get_actor(grid))

        #============== Boundary types
        if self.bnd_view.model() is not self.bnd_model:
            self.bnd_view.setModel(self.bnd_model)
        self.bnd_model.reset()

        #============== CONTOURS
        if self.cont_view.model() is not self.cont_model:
            self.cont_view.setModel(self.cont_model)
        self.cont_model.reset()
        self.cont_view.expandAll()

        #============== Draw
        #set actors list
        alist = grid_actors + self.cont_model.active_actors()
        self.vtk_draw.set_actor_list(alist)
        self.vtk_draw.ren.ResetCamera()
        self.vtk_draw.Render()
