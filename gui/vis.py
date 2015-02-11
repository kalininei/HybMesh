import sip

from PyQt4.QtGui import QTreeWidgetItem, QIcon, QPixmap
import vtk

import gridcom
import globvars
import dlgs

def_color = (0.7, 0.7, 0.7)
def_linewidth = 1
chk_color = (1, 0, 0)
chk_linewidth = 2


class FrameworkVis(object):
    'Interface for visual representation of framework object'
    def __init__(self):
        self.fw = None

    def set_framework(self, fw):
        ' sets a refference to the framework object @fw '
        self.fw = fw

    def update(self):
        ' Updates view of referenced framework '
        raise NotImplementedError


class GridBlockTreeItem(QTreeWidgetItem):
    ' presents a gridblock item in a grid manager docker window '
    __icons = {}

    @classmethod
    def get_icon(cls, code):
        ' get icon by its string code '
        if len(cls.__icons) == 0:
            cls.__icons['eye-on'] = QIcon(QPixmap(":/icons/eye-on.png"))
            cls.__icons['eye-off'] = QIcon(QPixmap(":/icons/eye-off.png"))
            cls.__icons['opts'] = QIcon(QPixmap(":/icons/opts.png"))
            cls.__icons['del'] = QIcon(QPixmap(":/icons/delete.png"))
        return cls.__icons[code]

    def __init__(self, name, view_opt):
        super(GridBlockTreeItem, self).__init__()
        self.view_opt = view_opt
        self.name = name
        #options
        self.setIcon(3, self.get_icon("opts"))
        #delete button
        self.setIcon(4, self.get_icon("del"))
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
            self.setIcon(2, self.get_icon("eye-on"))
        else:
            self.setIcon(2, self.get_icon("eye-off"))

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
            com = gridcom.RemoveGrid2(self.name)
            globvars.actual_flow().exec_command(com)


class ViewGridData(object):
    """ A container for grid block representation
    Holds actors, Grid Manager QTreeWidgetItem

    All objects are unique for each grid_id and
    don't free memory with the deletion of the corresponding grid
    in order to support undo/redo strategy.
    Hence corresponding actors and other memory consuming attrs
    should be freed mannualy with the deletion of the grid by
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
        g2.add_subscriber_change_geom(ret)
        return ret

    def mem_minimize(self):
        " frees all memory consuming data "
        self._grid_actor = None
        self._manager_item = None

    def grid_geom_event(self, g2):
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
            map(cell.GetPointIds().InsertNextId, cl)
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
    """ output into two Qt Frames:
         -- visual
         -- grid tree
    """
    def __init__(self):
        super(QtFrameworkVis, self).__init__()

    def get_checked_grid_names(self):
        ' returns names of checked grids '
        ret = []
        for (name, grid) in self.fw.grids2.items():
            view_data = ViewGridData(grid)
            if view_data.get_checked():
                ret.append(name)
        return ret

    def update(self):
        """ updates view in of referenced framework in the main window '
            after new command has been executed
        """
        # cold start invocation before window.show() call
        if not globvars.mainWindow.isVisible():
            return

        #view grids update
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

        #set actors list
        globvars.mainWindow.set_grid_actors(grid_actors)

        #draw
        globvars.mainWindow.vtkWidget.GetRenderWindow().Render()
        globvars.mainWindow.vtkWidget.ren.ResetCamera()
