#!/usr/bin/env python
'geometric picker dialogs'
import copy
from PyQt4 import QtGui, QtCore
import vtk
import globvars
import drawarea
import dlgs
import bgeom
import optview


class _SubPoly(object):
    """ vtkPolydata which is constructed from the subset of
        contours edges
    """
    def __init__(self, conts):
        '(Contour2)'
        self.conts = conts
        self.start_ed = [0]
        self.start_pnt = [0]
        for c in self.conts:
            self.start_ed.append(
                    self.start_ed[-1] + c.n_edges())
            self.start_pnt.append(
                    self.start_pnt[-1] + c.n_points())

        #build points
        ps = vtk.vtkPoints()
        ps.Allocate(self.start_pnt[-1])
        for c in self.conts:
            for p in c.points:
                ps.InsertNextPoint(p.x, p.y, 0)

        #build poly
        self.polydata = vtk.vtkPolyData()
        self.polydata.SetPoints(ps)

        mp = vtk.vtkPolyDataMapper()
        mp.SetInputData(self.polydata)
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(mp)

        #actor properties
        chk_color = [x / 256.0 for x in globvars.view_options.checked_color]
        chk_linewidth = globvars.view_options.checked_cont_line_width
        self.actor.GetProperty().SetColor(*chk_color)
        self.actor.GetProperty().SetLineWidth(chk_linewidth)

    def set_eds_list(self, edsets):
        '[set-of-edges]'
        self.polydata.DeleteCells()
        eds = vtk.vtkCellArray()
        for cind, s in enumerate(edsets):
            if len(s) == 0:
                continue
            ep = self.conts[cind].edges_points()
            for e in s:
                ed = vtk.vtkLine()
                ed.GetPointIds().SetId(0, self._glob_pnt(cind, ep[e][0]))
                ed.GetPointIds().SetId(1, self._glob_pnt(cind, ep[e][1]))
                eds.InsertNextCell(ed)
        self.polydata.SetLines(eds)

    def _glob_pnt(self, cont_ind, p_ind):
        return self.start_pnt[cont_ind] + p_ind


class _BndModel(QtCore.QAbstractTableModel):
    def __init__(self, bnds, parent=None):
        super(_BndModel, self).__init__(parent)
        self.bnds = bnds

    def columnCount(self, parent=None):
        return 2

    def rowCount(self, parent=None):
        return self.bnds.bnd_count() + 1

    def data(self, index, role):
        bind = index.row()
        if bind == self.rowCount() - 1:
            bind = -1
        col = index.column()
        n = self.bnds.get_names()[bind]
        b = self.bnds.get(name=n)
        if role == QtCore.Qt.DisplayRole:
            if bind < 0 and col == 0:
                return "Initial values"
            if col == 0 and bind >= 0:
                return "%i: %s" % (b.index, b.name)
        if role == QtCore.Qt.BackgroundRole:
            if bind >= 0 and col == 1:
                return QtGui.QColor(*b.color)
        if role == QtCore.Qt.SizeHintRole:
            if col == 1:
                return QtCore.QSize(20, 15)

    def flags(self, index):
        ret = QtCore.Qt.ItemIsEnabled
        if index.column() != 1:
            ret |= QtCore.Qt.ItemIsSelectable
        return ret


class _BndList(QtGui.QTableView):
    def __init__(self, bnds, parent=None):
        super(_BndList, self).__init__(parent)
        self.setModel(_BndModel(bnds, parent))
        hh = self.horizontalHeader()
        hh.setVisible(False)
        hh.setResizeMode(QtGui.QHeaderView.ResizeToContents)
        vh = self.verticalHeader()
        vh.setResizeMode(QtGui.QHeaderView.ResizeToContents)
        vh.setVisible(False)
        self.setShowGrid(False)
        self.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.selectRow(0)

        self.bnds = bnds

    def mousePressEvent(self, event):
        #don't lose selection on outside item press
        index = self.indexAt(event.pos())
        if index.isValid() and index.column() == 0:
            super(_BndList, self).mousePressEvent(event)

    def current_b(self):
        """ returns index of currently selected boundary
            returns -1 for initial state check
        """
        row = self.selectedIndexes()[0].row()
        if row == self.model().rowCount() - 1:
            ret = -1
        else:
            ret = self.bnds.get(orderindex=row).index
        return ret


class AskForRect(dlgs.SimpleAbstractDialog):
    def __init__(self, parent=None):
        super(AskForRect, self).__init__(parent)
        self.resize(300, 200)
        self.setWindowTitle("Enter selecion rectangle")

    def _default_odata(self, obj):
        "-> options struct with default values"
        obj.p0, obj.p1 = bgeom.Point2(0.0, 0.0), bgeom.Point2(1.0, 1.0)

    def olist(self):
        "-> optview.OptionsList"
        return optview.OptionsList([
            ("Rectangle", "Bottom left",
                optview.XYOptionEntry(self.odata(), "p0")),
            ("Rectangle", "Top right",
                optview.XYOptionEntry(self.odata(), "p1"))])

    def ret_value(self):
        ' -> p0, p1 '
        od = copy.deepcopy(self.odata())
        if od.p0.x > od.p1.x:
            od.p0.x, od.p1.x = od.p1.x, od.p0.x
        if od.p0.y > od.p1.y:
            od.p0.y, od.p1.y = od.p1.y, od.p0.y
        return od.p0, od.p1


class BoundaryTypesManagement(QtGui.QDialog):
    #selection behaviours types
    BEHAVIOUR_AEXCL = 0
    BEHAVIOUR_UNITE = 1
    BEHAVIOUR_SEPAR = 2

    def __init__(self, conts, bnd_types, parent=None):
        """
            conts = [(name, Contour2, vtkPolyData), ...]
            bnd_types = btypes.BndTypesList
        """
        #================ Representation
        super(BoundaryTypesManagement, self).__init__(parent)
        self.setStyleSheet("QGroupBox { font-weight: bold; } ")
        self.setWindowTitle("Assign boundary types to contours")
        self.setModal(True)
        #==== basic layout
        splitter = QtGui.QSplitter(QtCore.Qt.Horizontal, self)
        splitter.setChildrenCollapsible(False)
        layout = QtGui.QHBoxLayout(self)
        layout.addWidget(splitter)

        #==== left frame
        self.wdraw = drawarea.VTKWidget(splitter)

        #==== right frame
        self.inp = QtGui.QFrame(self)

        #set bc widgets
        self.frame_set = QtGui.QGroupBox("Set boundary type", self)
        self.tab_bc = _BndList(bnd_types, self)
        self.but_set_sel = QtGui.QPushButton("To selected")
        self.but_set_all = QtGui.QPushButton("To all")

        setlayout = QtGui.QGridLayout(self.frame_set)
        setlayout.addWidget(self.tab_bc, 0, 0, 1, 2)
        setlayout.addWidget(self.but_set_sel, 1, 0)
        setlayout.addWidget(self.but_set_all, 1, 1)

        #selection widget
        self.frame_select = QtGui.QGroupBox("Selection", self)
        self.rb_excluding = QtGui.QRadioButton("Auto excluding")
        self.rb_union = QtGui.QRadioButton("Union")
        self.rb_union.setChecked(True)
        self.rb_separated = QtGui.QRadioButton("Separated")
        self.but_reverse = QtGui.QPushButton("Reverse")
        self.but_clear = QtGui.QPushButton("Clear")
        self.but_rect = QtGui.QPushButton("Rect ...")

        sellayout = QtGui.QGridLayout(self.frame_select)
        sellayout.addWidget(self.but_clear, 0, 0)
        sellayout.addWidget(self.but_reverse, 0, 1)
        sellayout.addWidget(self.but_rect, 1, 0)
        sellayout.addWidget(QtGui.QLabel("Behaviour:"), 2, 0, 1, 2)
        sellayout.addWidget(self.rb_excluding, 3, 0, 1, 2)
        sellayout.addWidget(self.rb_union, 4, 0, 1, 2)
        sellayout.addWidget(self.rb_separated, 5, 0, 1, 2)

        #ok/cancel
        bts = QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel
        ori = QtCore.Qt.Horizontal
        bbox = QtGui.QDialogButtonBox(bts, ori)
        bbox.accepted.connect(self.accept)
        bbox.rejected.connect(self.reject)

        inplayout = QtGui.QVBoxLayout(self.inp)
        inplayout.addWidget(self.frame_set)
        inplayout.addWidget(self.frame_select)
        inplayout.addWidget(bbox)

        #==== main positioning
        splitter.addWidget(self.wdraw)
        splitter.addWidget(self.inp)
        splitter.setStretchFactor(1, 0)
        splitter.setStretchFactor(0, 1)
        splitter.setSizes([0, 200])

        self._assign_buttons()

        #================== Init
        self.bt = bnd_types
        self.names = [x[0] for x in conts]
        self.conts = [x[1] for x in conts]
        self.polys = [vtk.vtkPolyData() for _ in conts]
        for p, x in zip(self.polys, conts):
            p.DeepCopy(x[2])
        self.old_bnd = [[] for _ in conts]
        self.new_bnd = [[] for _ in conts]
        self.col_bnd = [vtk.vtkUnsignedCharArray() for _ in conts]
        #[set-of-edge-indicies]
        self.selected_cells = [set() for _ in conts]

        self.checked_poly = _SubPoly(self.conts)
        self.acts = []
        for i in range(len(self.names)):
            self.acts.append(self._init_data(i))
        self.wdraw.set_actor_list(self.acts + [self.checked_poly.actor])
        self._start_selection_ineractor()

    def _assign_buttons(self):
        'buttons events assign'
        self.but_set_sel.clicked.connect(self._set_to_selected)
        self.but_set_all.clicked.connect(self._set_all_to_selected)
        self.but_reverse.clicked.connect(self._reverse_selection)
        self.but_clear.clicked.connect(self._clear_all_selection)
        self.but_rect.clicked.connect(self._select_rect_dialog)

    def _get_behaviour(self):
        'returns selection behaviour'
        if self.rb_excluding.isChecked():
            return self.BEHAVIOUR_AEXCL
        if self.rb_union.isChecked():
            return self.BEHAVIOUR_UNITE
        if self.rb_separated.isChecked():
            return self.BEHAVIOUR_SEPAR

    def _set_to_selected(self):
        'assign selected edges with self.selected_cells'
        b = self.tab_bc.current_b()
        for i, st in enumerate(self.selected_cells):
            if len(st) > 0:
                self._set_boundary(i, b, st)
        self._clear_all_selection()

    def _set_all_to_selected(self):
        for ind, c in enumerate(self.conts):
            self.selected_cells[ind].update(range(c.n_edges()))
        self._set_to_selected()

    def _reverse_selection(self):
        for ind, c in enumerate(self.conts):
            self._add_selection(ind, range(c.n_edges()), self.BEHAVIOUR_AEXCL)

    def _select_rect_dialog(self):
        dialog = AskForRect()
        if dialog.exec_():
            ret = dialog.ret_value()
            self._select_in_rectangle(ret[0].x, ret[0].y, ret[1].x, ret[1].y)

    def _select_in_rectangle(self, xmin, ymin, xmax, ymax):
        'select area in world coordinates'
        dc1, dc2 = [0, 0, 0], [0, 0, 0]
        vtk.vtkInteractorObserver.ComputeWorldToDisplay(
                self.wdraw.ren, xmin, ymin, 0, dc1)
        vtk.vtkInteractorObserver.ComputeWorldToDisplay(
                self.wdraw.ren, xmax, ymax, 0, dc2)
        self._select_in_display_rect(dc1[0], dc1[1], dc2[0], dc2[1])

    def _select_in_display_rect(self, xmin, ymin, xmax, ymax):
        'select area in display coordinates'
        #removes checked actor to exclude it from selection
        self.wdraw.remove_actor(self.checked_poly.actor)

        selector = vtk.vtkHardwareSelector()
        selector.SetFieldAssociation(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS)
        selector.SetRenderer(self.wdraw.ren)

        selector.SetArea(int(xmin), int(ymin), int(xmax), int(ymax))
        selection = selector.Select()

        num = selection.GetNumberOfNodes()
        s_edges = [[] for _ in self.conts]
        if num > 0:
            for i in range(num):
                sel_node = selection.GetNode(i)
                dt = sel_node.GetSelectionList()
                n = dt.GetNumberOfTuples()
                #indicate selected contour index
                actor = sel_node.GetProperties().Get(
                        vtk.vtkSelectionNode.PROP())
                cindex = self.acts.index(actor)
                s_edges[cindex] = [dt.GetValue(i) for i in range(n)]
        for i, s in enumerate(s_edges):
            self._add_selection(i, s)

        #add selection actor back
        self.wdraw.add_actor(self.checked_poly.actor)

    def _init_data(self, ind):
        '->vtkActor'
        self.new_bnd[ind] = [0] * self.conts[ind].n_edges()
        self.col_bnd[ind].SetNumberOfComponents(3)
        self.col_bnd[ind].SetNumberOfValues(3 * self.conts[ind].n_edges())
        self.col_bnd[ind].SetName("__temp_col")
        for i in range(self.conts[ind].n_edges()):
            b = self.conts[ind].edge_bnd(i)
            self.new_bnd[ind][i] = b
        self.old_bnd = copy.deepcopy(self.new_bnd)
        self.polys[ind].GetCellData().AddArray(self.col_bnd[ind])
        self._update_colors(ind)

        mp = vtk.vtkPolyDataMapper()
        mp.SetInputData(self.polys[ind])
        actor = vtk.vtkActor()
        actor.SetMapper(mp)
        return actor

    def _update_colors(self, ind):
        'updates colors according to self.new_bnd array'
        for i, b in enumerate(self.new_bnd[ind]):
            self.col_bnd[ind].SetTupleValue(i, self.bt.get(index=b).color)
        self.polys[ind].GetCellData().SetActiveScalars("__temp_col")
        self.wdraw.Render()

    def _set_boundary(self, ind, bnd_index, edges):
        'sets bnd_index boundary to edges of conts[ind]'
        for e in edges:
            if bnd_index < 0:
                b = self.old_bnd[ind][e]
            else:
                b = bnd_index
            self.new_bnd[ind][e] = b
        self._update_colors(ind)

    def _clear_all_selection(self):
        for s in self.selected_cells:
            s.clear()
        #draw selection
        self.checked_poly.set_eds_list(self.selected_cells)
        self.wdraw.Render()

    def _add_selection(self, ind, edges, behav=None):
        if behav is None:
            behav = self._get_behaviour()
        s = self.selected_cells[ind]
        if behav == self.BEHAVIOUR_AEXCL:
            s.symmetric_difference_update(edges)
        elif behav == self.BEHAVIOUR_SEPAR:
            s.clear()
            s.update(edges)
        elif behav == self.BEHAVIOUR_UNITE:
            s.update(edges)
        #draw selection
        self.checked_poly.set_eds_list(self.selected_cells)
        self.wdraw.Render()

    def exec_(self):
        'overriden'
        self.show()
        self.wdraw.Initialize()
        return super(BoundaryTypesManagement, self).exec_()

    def ret_value(self):
        '->[{boundary-index: [list-of-edges] }, ... same for all contours]'
        ret = []
        for i, c in enumerate(self.conts):
            ret.append({})
            r = ret[-1]
            for s in set(self.new_bnd[i]):
                r[s] = []
            for i, e in enumerate(self.new_bnd[i]):
                r[e].append(i)
        return ret

    def _start_selection_ineractor(self):
        inter = vtk.vtkInteractorStyleRubberBand2D()
        inter.AddObserver(vtk.vtkCommand.SelectionChangedEvent,
                self._mouse_selection)
        self.wdraw.iren().SetInteractorStyle(inter)

    def _mouse_selection(self, caller, event):
        'callback after mouse selection'
        xmin, ymin = caller.GetStartPosition()
        xmax, ymax = caller.GetEndPosition()
        if xmin > xmax:
            xmin, xmax = xmax, xmin
        if ymin > ymax:
            ymin, ymax = ymax, ymin
        self._select_in_display_rect(xmin, ymin, xmax, ymax)


if __name__ == "__main__":
    pass
