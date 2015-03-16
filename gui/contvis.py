"contours visualization"

from PyQt4 import QtCore, QtGui
import vtk
import qtbp
import globvars
import objcom


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
        if not index.isValid():
            return
        if event.button() == QtCore.Qt.RightButton:
            pass
        if event.button() == QtCore.Qt.LeftButton:
            #remove user contour
            if index.column() == 3 and \
                    index.parent().data(0) == "User defined contours":
                com = objcom.RemoveGrid([], [index.internalPointer().data(0)])
                globvars.actual_flow().exec_command(com)
            #switches
            else:
                self.model().switch_item_status(index)


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
            #delete exciting actor on contour contour geometry change
            ret.contour.add_subscriber_change_geom(ret._delete_actor)
            #actor
            ret._actor = None
            #view options
            ret._ischecked = False
            ret._isvisible = True
            #backup created object
            cls.stored[i] = ret
        ret.name = name
        return ret

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

    def switch_check(self):
        'switches check status of the contour'
        self._ischecked = not self._ischecked
        self._set_view()

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

    def switch_visibility(self):
        'switches visibility status of the contour'
        self._isvisible = not self._isvisible
        if self._isvisible:
            self.get_actor().VisibilityOn()
        else:
            self.get_actor().VisibilityOff()


class ContListModel(qtbp.RowsTreeModel):
    """ Model for contours list representation
    """

    def __init__(self, clist, glist, b, parent=None):
        """ ({name -> contours2.Contour} clist,
             {name -> grid2.Grid2} glist,
             btypes.BndTypesList b, QWidget parent)"
        """
        super(ContListModel, self).__init__(UserContourEntry.col_count, parent)
        self.uclist = clist
        self.glist = glist
        self.bnd = b

    def contitem_by_index(self, index):
        'returns existing UserContourEntry from index'
        item = index.internalPointer()
        if item.parent() is self.user_cont:
            cont = self.uclist[item.data(0)]
        elif item.parent() is self.grid_cont:
            cont = self.glist[item.data(0)].cont
        else:
            return None
        return UserContourEntry(None, cont)

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

    def get_actor(self, name):
        '-> vtkActor. Get actor by contour name'
        try:
            cont = self.uclist[name]
        except KeyError:
            cont = self.glist[name].cont
        v = UserContourEntry(None, cont)
        return v.get_actor()

    def reset(self):
        "overriden"
        self.clist = [(k, v) for k, v in self.uclist.items()]
        self.clist += [(k, v.cont) for k, v in self.glist.items()]

        #remove all data from tree representation
        #and rebuild data
        self._root_item.clear_childs()
        self._fill_data_tree(self.uclist, self.glist)
        #set colors to actors
        for v in self.clist:
            UserContourEntry(*v).set_boundary_type(self.bnd)
        super(ContListModel, self).reset()

    def data(self, index, role):
        "overriden"
        item = index.internalPointer()
        #rows captions
        if index.column() == 0 and role == QtCore.Qt.DisplayRole:
            return item.data(0)
        #buttons
        if item.parent() is not self._root_item:
            cont_item = self.contitem_by_index(index)
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

    def switch_item_status(self, index):
        "this is called when left mouse button event occurs"
        item = index.internalPointer()
        if item.parent() is self._root_item:
            return
        if index.column() == 0:
            return
        cont_item = self.contitem_by_index(index)
        #switch check status
        if index.column() == 1:
            cont_item.switch_check()
        if index.column() == 2:
            cont_item.switch_visibility()
        self.dataChanged.emit(index, index)
        globvars.mainvtk_render()

    def active_actors(self):
        '-> list of vtk actors'
        return [UserContourEntry(*v).get_actor() for v in self.clist]

    def _fill_data_tree(self, clist, glist):
        #grid contours
        self.grid_cont = qtbp.TreeItem(
                self._treeitem_data("Grids boundaries"), self._root_item)
        for k, v in glist.items():
            e = UserContourEntry(k, v.cont)
            qtbp.TreeItem(e.full_data(), self.grid_cont)

        #user contours
        self.user_cont = qtbp.TreeItem(
                self._treeitem_data("User defined contours"), self._root_item)
        for k, v in clist.items():
            e = UserContourEntry(k, v)
            qtbp.TreeItem(e.full_data(), self.user_cont)
