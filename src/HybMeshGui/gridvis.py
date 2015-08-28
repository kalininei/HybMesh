"grids visualization"

from PyQt4 import QtCore, QtGui
import vtk
import qtbp
import globvars
import objcom
import functools
import gridcom


class GridManagerItemView(QtGui.QTreeView):
    def __init__(self, parent=None):
        super(GridManagerItemView, self).__init__(parent)
        self.setHeaderHidden(True)
        self.setSelectionBehavior(QtGui.QAbstractItemView.SelectItems)
        self.header().setResizeMode(
                QtGui.QHeaderView.ResizeToContents)
        self.header().setStretchLastSection(False)

    def mousePressEvent(self, event):
        'overriden'
        index = self.indexAt(event.pos())
        if event.button() == QtCore.Qt.RightButton:
            #context menu
            acts = self.model().context_menu_items(index, self)
            if acts is not None:
                menu = QtGui.QMenu(self)
                map(menu.addAction, acts)
                menu.popup(self.viewport().mapToGlobal(event.pos()))

        if event.button() == QtCore.Qt.LeftButton:
            #left mouse click
            self.model().mouse_click(index)


class GridEntry(object):
    """ grid representation:
        vtk actor + visibility options
    """
    #geom id -> GridEntry
    stored = {}
    #number of column for list representation
    col_count = 5

    def __new__(cls, name, grid):
        i = grid.get_id()
        if i in cls.stored:
            ret = cls.stored[i]
        else:
            ret = super(GridEntry, cls).__new__(cls)
            ret.grid = grid
            #delete exciting actor on contour contour geometry change
            ret.grid.add_subscriber_change_geom(ret._delete_actor)
            #actor
            ret._actor = None
            #view options
            ret._ischecked = False
            ret._isvisible = True
            #backup created object
            cls.stored[i] = ret
        if name is not None:
            ret.name = name
        return ret

    @classmethod
    def mem_minimize(cls, actual_data):
        """ clears all stored geometry data for
            all objects except actual_data
        """
        for v in cls.stored.values():
            if v.grid not in actual_data:
                v._delete_actor()

    def full_data(self):
        'return data as a tuple for qtbp.TreeItem constructor'
        return (self.name, None, None, None, None)

    def get_actor(self):
        '->vtkActor'
        if self._actor is None:
            self._reset_actor()
            self._set_view()
        return self._actor

    def _delete_actor(self):
        self._actor = None
        self._poly = None

    def _reset_actor(self):
        """ Create an actor for the contour.
            Fills self._actor and self._poly fields.
        """
        ps = vtk.vtkPoints()
        poly = vtk.vtkUnstructuredGrid()

        #points
        ps.SetNumberOfPoints(self.grid.n_points())
        for i, pnt in enumerate(self.grid.points):
            ps.InsertPoint(i, pnt.x, pnt.y, 0)
        poly.SetPoints(ps)

        #cells
        for i, cl in enumerate(self.grid.cells_nodes_connect()):
            cell = vtk.vtkPolygon()
            cell.GetPointIds().Allocate(len(cl))
            #WTF? without reversing some edges are black
            map(cell.GetPointIds().InsertNextId, cl[::-1])
            poly.InsertNextCell(cell.GetCellType(), cell.GetPointIds())

        #mapper and actor
        mp = vtk.vtkDataSetMapper()
        mp.SetInputData(poly)
        self._actor = vtk.vtkActor()
        self._actor.GetProperty().SetRepresentationToWireframe()
        self._actor.SetMapper(mp)
        self._poly = poly

    def _set_view(self):
        'sets self._actor line width and color'
        if self._ischecked:
            w = globvars.view_options.checked_grid_line_width
            color = globvars.view_options.checked_color
        else:
            w = globvars.view_options.grid_line_width
            color = globvars.view_options.grid_color
        self._actor.GetProperty().SetLineWidth(w)
        color = tuple(c / 255.0 for c in color)
        self._actor.GetProperty().SetColor(*color)
        if self._isvisible:
            self.get_actor().VisibilityOn()
        else:
            self.get_actor().VisibilityOff()

    def switch_visibility(self, force=None):
        'switches visibility status of the contour. force - True/False'
        if force is not None and self._isvisible == force:
            return
        self._isvisible = not self._isvisible
        self._set_view()

    def switch_check(self, force=None):
        'switches check status of the contour'
        if force is not None and self._ischecked == force:
            return
        self._ischecked = not self._ischecked
        self._set_view()


class GridListModel(qtbp.RowsTreeModel):
    """ Model for grid list representation
    """

    def __init__(self, fwvis, parent=None):
        """ (vis.QtFrameworkVis fwvis, QWidget parent)
        """
        super(GridListModel, self).__init__(GridEntry.col_count, parent)
        self.fwvis = fwvis
        self.glist = fwvis.fw.grids2

    def _fill_data_tree(self, glist):
        #grids
        for k, v in glist.items():
            e = GridEntry(k, v)
            qtbp.TreeItem(e.full_data(), self._root_item)

    # ================ easy data access
    def _griditem_by_index(self, index):
        'returns existing GridEntry from Qt index'
        item = index.internalPointer()
        g = self.glist[item.data(0)]
        return GridEntry(None, g)

    def _all_entries(self):
        '->[GridEntry] for all contours'
        return [GridEntry(None, g) for g in self.glist.values()]

    # ================= ItemModel behavior
    def reset(self):
        "overriden"
        #clears geometry for all except active data
        GridEntry.mem_minimize(self.glist.values())

        #remove all data from tree representation
        #and rebuild data
        self._root_item.clear_childs()
        self._fill_data_tree(self.glist)
        super(GridListModel, self).reset()

    def data(self, index, role):
        "overriden"
        item = index.internalPointer()
        #rows captions
        if index.column() == 0 and role == QtCore.Qt.DisplayRole:
            return item.data(0)
        #buttons
        grid_item = self._griditem_by_index(index)
        #checkbox
        if index.column() == 1 and role == QtCore.Qt.CheckStateRole:
            return 2 if grid_item._ischecked else 0
        #visible
        if index.column() == 2 and role == QtCore.Qt.DecorationRole:
            if grid_item._isvisible:
                return qtbp.get_icon("eye-on")
            else:
                return qtbp.get_icon("eye-off")
        #options
        if index.column() == 3 and role == QtCore.Qt.DecorationRole:
            return qtbp.get_icon("opts")
        #remove
        if index.column() == 4 and role == QtCore.Qt.DecorationRole:
            return qtbp.get_icon("del")
        #size
        if role == QtCore.Qt.SizeHintRole:
            if index.column() > 0:
                return QtCore.QSize(20, 20)

    def flags(self, index):
        "overriden"
        ret = QtCore.Qt.ItemIsEnabled
        return ret

    def _reset_table(self):
        'resets icons in the view table'
        ind1 = self.index(0, 0, QtCore.QModelIndex())
        ind2 = self.index(len(self.glist) - 1, GridEntry.col_count - 1,
                QtCore.QModelIndex())
        self.dataChanged.emit(ind1, ind2)

    # ======================== Mouse Interactions
    def mouse_click(self, index):
        "this is called when left mouse button event occurs"
        if not index.isValid():
            return
        item = index.internalPointer()
        if index.column() == 0:
            return
        elif index.column() == 4:
            #remove user contour
            com = objcom.RemoveGrid([item.data(0)], [])
            globvars.actual_flow().exec_command(com)
        elif index.column() == 3:
            #view opt
            pass
        elif index.column() in [1, 2]:
            grid_item = self._griditem_by_index(index)
            if index.column() == 2:
                #switch visibility
                grid_item.switch_visibility()
            elif index.column() == 1:
                #switch check
                grid_item.switch_check()
            self.dataChanged.emit(index, index)
            globvars.mainvtk_render()

    def _show_all(self, yesno):
        [e.switch_visibility(yesno) for e in self._all_entries()]
        globvars.mainvtk_render()
        self._reset_table()

    def _check_all(self, yesno):
        [e.switch_check(yesno) for e in self._all_entries()]
        globvars.mainvtk_render()
        self._reset_table()

    def _show_only(self, index):
        nm = str(index.data().toString())
        for e in self._all_entries():
            e.switch_visibility(e.name == nm)
        self.fwvis.show_cont(None, False)

    def _check_only(self, index):
        nm = str(index.data().toString())
        for e in self._all_entries():
            e.switch_check(e.name == nm)
        self.fwvis.check_cont(None, False)

    def _show_grid_and_contour(self, index):
        nm = str(index.data().toString())
        self.fwvis.show_cont(None, False)
        self.fwvis.show_grid(None, False)
        self.fwvis.show_cont(nm, True)
        self.fwvis.show_grid(nm, True)

    def _rename(self, index):
        oldname = index.data().toString()
        newname, ok = QtGui.QInputDialog.getText(None,
                "Rename grid", "Enter new grid name", text=oldname)
        oldname, newname = str(oldname), str(newname)
        if ok and oldname != newname:
            com = gridcom.RenameGrid(oldname, newname)
            globvars.actual_flow().exec_command(com)

    def context_menu_items(self, index, parent):
        """ (QModelIndex, QObject) -> [QtGui.QAction].
            Gives a context menu for the item by its Qt index.
        """
        ret = []

        def add_act(name, func=None, *args):
            a = QtGui.QAction(name, parent)
            if func is not None:
                a.triggered.connect(functools.partial(func, *args))
            else:
                a.setSeparator(True)
            ret.append(a)

        if index.isValid():
            if index.column() != 0:
                return
            add_act("Show only", self._show_only, index)
            add_act("Check only", self._check_only, index)
            add_act("Show grid && cont", self._show_grid_and_contour, index)
            add_act("")
            add_act("Rename", self._rename, index)
        else:
            add_act("Show All", self._show_all, True)
            add_act("Hide All", self._show_all, False)
            add_act("")
            add_act("Check All", self._check_all, True)
            add_act("Uncheck All", self._check_all, False)
        return ret

    # =============== Methods for outer world
    def checked_names(self):
        '->[list-of-str] names of all checked grids'
        ret = []
        for k, v in self.glist.items():
            itm = GridEntry(None, v)
            if itm._ischecked:
                ret.append(k)
        return ret

    def check(self, grid, yesno):
        'Checks/Unchecks the grid2 object entry'
        if grid is not None:
            GridEntry(None, grid).switch_check(yesno)
        else:
            [e.switch_check(yesno) for e in self._all_entries()]
        self._reset_table()
        globvars.mainvtk_render()

    def show(self, grid, yesno):
        'Shows/Hides the grid2 object entry'
        if grid is not None:
            GridEntry(None, grid).switch_visibility(yesno)
        else:
            [e.switch_visibility(yesno) for e in self._all_entries()]
        self._reset_table()
        globvars.mainvtk_render()

    def active_actors(self):
        '-> list of vtk actors'
        return [GridEntry(*v).get_actor() for v in self.glist.items()]
