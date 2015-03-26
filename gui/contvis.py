"contours visualization"

from PyQt4 import QtCore, QtGui
import vtk
import qtbp
import globvars
import objcom
import contcom
import functools
import dlgs


class ContManagerItemView(QtGui.QTreeView):
    def __init__(self, parent=None):
        super(ContManagerItemView, self).__init__(parent)
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


class UserContourEntry(object):
    """ contour representation:
        vtk actor + visibility options
    """
    #cont id -> UserContourEntry
    stored = {}
    #number of column for list representation
    col_count = 4

    def __new__(cls, name, contour):
        i = contour.get_id()
        if i in cls.stored:
            ret = cls.stored[i]
        else:
            ret = super(UserContourEntry, cls).__new__(cls)
            ret.contour = contour
            #delete existing actor on contour contour geometry change
            ret.contour.add_subscriber_change_geom(ret._delete_actor)
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
            if v.contour not in actual_data:
                v._delete_actor()

    def set_boundary_type(self, bt):
        'sets a boundary type to contour'
        self.bt = bt
        self.get_actor()
        _actor_colors = vtk.vtkUnsignedCharArray()
        _actor_colors.SetNumberOfComponents(3)
        _actor_colors.SetName("Colors")
        _actor_colors.SetNumberOfValues(3 * self.contour.n_edges())
        for i in range(self.contour.n_edges()):
            b = self.contour.edge_bnd(i)
            _actor_colors.SetTupleValue(i, self.bt.get(index=b).color)
        self._poly.GetCellData().AddArray(_actor_colors)
        self._set_view()

    def full_data(self):
        'return data as a tuple for qtbp.TreeItem constructor'
        return (self.name, None, None, None)

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
        polydata = vtk.vtkPolyData()

        #points
        ps.SetNumberOfPoints(self.contour.n_points())
        for i, pnt in enumerate(self.contour.points):
            ps.InsertPoint(i, pnt.x, pnt.y, 0)
        polydata.SetPoints(ps)

        #cells
        eds = vtk.vtkCellArray()
        for i, cl in enumerate(self.contour.edges_points()):
            ed = vtk.vtkLine()
            ed.GetPointIds().SetId(0, cl[0])
            ed.GetPointIds().SetId(1, cl[1])
            eds.InsertNextCell(ed)
        polydata.SetLines(eds)

        #cell colors allocation
        _checked_colors = vtk.vtkUnsignedCharArray()
        _checked_colors.SetNumberOfComponents(3)
        _checked_colors.SetName("CheckedColors")
        _checked_colors.SetNumberOfValues(3 * self.contour.n_edges())
        cc = globvars.view_options.checked_color
        for i in range(self.contour.n_edges()):
            _checked_colors.SetTupleValue(i, cc)
        polydata.GetCellData().AddArray(_checked_colors)

        #mapper and actor
        mp = vtk.vtkPolyDataMapper()
        mp.SetInputData(polydata)
        self._actor = vtk.vtkActor()
        self._actor.SetMapper(mp)
        self._poly = polydata

    def _set_view(self):
        'sets self._actor line width and color'
        if self._ischecked:
            w = globvars.view_options.checked_cont_line_width
            color = "CheckedColors"
        else:
            w = globvars.view_options.cont_line_width
            color = "Colors"
        self._actor.GetProperty().SetLineWidth(w)
        self._poly.GetCellData().SetActiveScalars(color)
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


class ContListModel(qtbp.RowsTreeModel):
    """ Model for contours list representation
    """

    def __init__(self, fwvis, parent=None):
        """ (vis.QtFrameworkVis fwvis, QWidget parent)
        """
        super(ContListModel, self).__init__(UserContourEntry.col_count, parent)
        self.fwvis = fwvis
        self.uclist = fwvis.fw.contours2
        self.glist = fwvis.fw.grids2
        self.bnd = fwvis.fw.boundary_types

    def _fill_data_tree(self, clist, glist):
        #grid contours
        self.grid_cont = qtbp.TreeItem(
                self._treeitem_data("Grids boundaries"), self._root_item)
        for k, v in glist.items():
            e = UserContourEntry(k, v.cont)
            qtbp.TreeItem(e.full_data(), self.grid_cont)

        #user contours
        self.user_cont = qtbp.TreeItem(
                self._treeitem_data("User contours"), self._root_item)
        for k, v in clist.items():
            e = UserContourEntry(k, v)
            qtbp.TreeItem(e.full_data(), self.user_cont)

    # ================== easy data access
    def _contitem_by_index(self, index):
        'returns existing UserContourEntry from Qt index'
        item = index.internalPointer()
        if item.parent() is self.user_cont:
            cont = self.uclist[item.data(0)]
        elif item.parent() is self.grid_cont:
            cont = self.glist[item.data(0)].cont
        else:
            return None
        return UserContourEntry(None, cont)

    def _all_entries(self):
        '->[UserContourEntry] for all contours'
        return [UserContourEntry(None, c) for c in self.uclist.values()] + \
            [UserContourEntry(None, g.cont) for g in self.glist.values()]

    def _cont_by_name(self, nm):
        '-> AbstractContour by string name'
        if nm in self.uclist:
            return self.uclist[nm]
        else:
            return self.glist[nm].cont

    # ============================ ItemModel behavior
    def reset(self):
        "overriden"
        self.clist = [(k, v) for k, v in self.uclist.items()]
        self.clist += [(k, v.cont) for k, v in self.glist.items()]

        #clears geometry for all except active data
        UserContourEntry.mem_minimize([c[1] for c in self.clist])

        #remove all data from tree representation
        #and rebuild data
        self._root_item.clear_childs()
        self._fill_data_tree(self.uclist, self.glist)
        #set colors to actors
        for v in self.clist:
            UserContourEntry(*v).set_boundary_type(self.bnd)
        super(ContListModel, self).reset()

    def _reset_table(self):
        'resets icons in the view table'
        ind1 = self.index(0, 0, QtCore.QModelIndex())
        ind2 = self.index(1, 0, QtCore.QModelIndex())
        self.dataChanged.emit(ind1, ind2)

    def data(self, index, role):
        "overriden"
        item = index.internalPointer()
        #rows captions
        if index.column() == 0 and role == QtCore.Qt.DisplayRole:
            return item.data(0)
        #buttons
        if item.parent() is not self._root_item:
            cont_item = self._contitem_by_index(index)
            #checkbox
            if index.column() == 1 and role == QtCore.Qt.CheckStateRole:
                return 2 if cont_item._ischecked else 0
            #visible
            if index.column() == 2 and role == QtCore.Qt.DecorationRole:
                if cont_item._isvisible:
                    return qtbp.get_icon("eye-on")
                else:
                    return qtbp.get_icon("eye-off")
            #remove
            if index.column() == 3 and item.parent() is self.user_cont \
                    and role == QtCore.Qt.DecorationRole:
                return qtbp.get_icon("del")
        #size
        if role in [QtCore.Qt.SizeHintRole]:
            if index.column() == 1:
                return QtCore.QSize(20, 20)
        #color
        if role in [QtCore.Qt.FontRole]:
            if item.parent() is self._root_item:
                f = QtGui.QFont()
                f.setBold(True)
                return f

    def flags(self, index):
        "overriden"
        ret = QtCore.Qt.ItemIsEnabled
        return ret

    # ======================== Mouse interaction
    def mouse_click(self, index):
        "this is called when left mouse button event occurs"
        if not index.isValid():
            return
        item = index.internalPointer()
        if item.parent() is self._root_item or index.column() == 0:
            return

        if index.column() == 3 and item.parent() == self.user_cont:
            #remove user contour
            com = objcom.RemoveGrid([], [item.data(0)])
            globvars.actual_flow().exec_command(com)
        else:
            cont_item = self._contitem_by_index(index)
            #switch check status
            if index.column() == 1:
                cont_item.switch_check()
            #switch visibility
            if index.column() == 2:
                cont_item.switch_visibility()
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
        self.fwvis.show_grid(None, False)

    def _check_only(self, index):
        nm = str(index.data().toString())
        for e in self._all_entries():
            e.switch_check(e.name == nm)
        self.fwvis.check_grid(None, False)

    def _show_grid_and_contour(self, index):
        nm = str(index.data().toString())
        self.fwvis.show_cont(None, False)
        self.fwvis.show_grid(None, False)
        self.fwvis.show_cont(nm, True)
        self.fwvis.show_grid(nm, True)

    def _set_bnd(self, index):
        item = self._contitem_by_index(index)
        ret = globvars.actual_data().ask_for_new_contours_bnd([item.name])
        if ret is not None:
            arg = [contcom.BTypePicker(item.name, ret[0])]
            com = contcom.SetBTypeToContour(arg)
            globvars.actual_flow().exec_command(com)

    def _simpl_cont(self, index):
        dialog = dlgs.SimplifyContour([str(index.data().toString())],
                self.uclist.keys())
        if dialog.exec_():
            v = dialog.ret_value()
            com = contcom.SimplifyContours(*v)
            globvars.actual_flow().exec_command(com)

    def _copy_to_user(self, index):
        nm = str(index.data().toString())
        all_grds = self.glist.keys()
        dialog = dlgs.GridBndToContour(nm, all_grds)
        if dialog.exec_():
            arg = dialog.ret_value()
            com = contcom.GridBndToContour(*arg)
            globvars.actual_flow().exec_command(com)

    def _rename(self, index):
        oldname = index.data().toString()
        newname, ok = QtGui.QInputDialog.getText(None,
                "Rename contour", "Enter new contour name", text=oldname)
        oldname, newname = str(oldname), str(newname)
        if ok and oldname != newname:
            com = contcom.RenameContour(oldname, newname)
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
            item = index.internalPointer()
            if index.column() != 0 or item.parent() is self._root_item:
                return
            add_act("Show only", self._show_only, index)
            add_act("Check only", self._check_only, index)
            if item.parent() is self.user_cont and index.column() == 0:
                add_act("")
                add_act("Rename", self._rename, index)
                add_act("Separate/simplify", self._simpl_cont, index)
            if item.parent() is self.grid_cont and index.column() == 0:
                add_act("Show grid && cont",
                        self._show_grid_and_contour, index)
                add_act("")
                add_act("Copy to user", self._copy_to_user, index)
            add_act("Set boundary type", self._set_bnd, index)
        else:
            add_act("Show All", self._show_all, True)
            add_act("Hide All", self._show_all, False)
            add_act("")
            add_act("Check All", self._check_all, True)
            add_act("Uncheck All", self._check_all, False)
        return ret

    # -------------------- Methods for outer world
    def checked_names(self):
        '->[list-of-str] names of all checked contours'
        return self.checked_names_udc() + self.checked_names_gbc()

    def checked_names_udc(self):
        '->[list-of-str] names of all checked user defined contours'
        ret = []
        for k, v in self.uclist.items():
            itm = UserContourEntry(None, v)
            if itm._ischecked:
                ret.append(k)
        return ret

    def checked_names_gbc(self):
        '->[list-of-str] names of all checked grid bnd contours'
        ret = []
        for k, v in self.glist.items():
            itm = UserContourEntry(None, v.cont)
            if itm._ischecked:
                ret.append(k)
        return ret

    def check(self, contour, yesno):
        'Checks/Unchecks the contour object entry'
        if contour is not None:
            UserContourEntry(None, contour).switch_check(yesno)
        else:
            [e.switch_check(yesno) for e in self._all_entries()]
        self._reset_table()
        globvars.mainvtk_render()

    def show(self, contour, yesno):
        'Shows/Hides the contour object entry'
        if contour is not None:
            UserContourEntry(None, contour).switch_visibility(yesno)
        else:
            [e.switch_visibility(yesno) for e in self._all_entries()]
        self._reset_table()
        globvars.mainvtk_render()

    def get_actor(self, contname):
        '-> vtkActor. Return vtkActor by contour name'
        return UserContourEntry(None, self._cont_by_name(contname)).get_actor()

    def active_actors(self):
        '-> list of vtk actors'
        return [UserContourEntry(*v).get_actor() for v in self.clist]
