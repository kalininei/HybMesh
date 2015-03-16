#!/usr/bin/env python

from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import (QMainWindow,
    QFileDialog, QHBoxLayout)
import qtui.ui_ComGridMain

import drawarea
import commandwin
import gridcom
import objcom
import contcom
import dlgs
import globvars
import contvis
import btypes


class MainWindow(QMainWindow, qtui.ui_ComGridMain.Ui_MainWindow):
    "Main program window"
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)

        self._connect_actions()
        self._grid_manager_init()
        self._cont_manager_init()
        self._bnd_manager_init()
        self._vtk_widget_init()

    def initialize(self):
        """ build window widgets:
            should be called after globvars.Flows init"""
        self._commands_history_init()
        globvars.Flows.add_actflow_subscriber(self)
        self.actual_flow_changed(globvars.Flows.get_actual_flow_name())

    # --------- constructor procedures
    def _connect_actions(self):
        "connects actions for main window"
        #file
        self.act_open.triggered.connect(self._open_flows)
        self.act_save.triggered.connect(self._save_flows)
        #edit
        self.act_undo.triggered.connect(self._undo_command)
        self.act_redo.triggered.connect(self._redo_command)
        #geometry
        self.act_unf_rect.triggered.connect(self._add_rectangle_grid)
        self.act_unf_circ.triggered.connect(self._add_circular_grid)
        self.act_unf_ring.triggered.connect(self._add_ring_grid)
        self.act_cont_rect.triggered.connect(self._add_rectangle_cont)
        self.act_movrot.triggered.connect(self._mov_rot)
        self.act_scale.triggered.connect(self._scale_grid)
        self.act_remove_grids.triggered.connect(self._remove_grids)
        self.act_copy_grids.triggered.connect(self._copy_grids)
        #tools
        self.act_unite_grids.triggered.connect(self._unite_grids)
        self.act_set_bc.triggered.connect(self._set_bc)

    def _grid_manager_init(self):
        ' initialization of the grid manager '
        self.tw_gridblocks = QtGui.QTreeWidget(self.DWGridManager)
        self.tw_gridblocks.resize(200, 500)
        self.tw_gridblocks.setSelectionMode(
                QtGui.QAbstractItemView.NoSelection)
        self.tw_gridblocks.setHeaderHidden(True)
        self.DWGridManager.setWidget(self.tw_gridblocks)
        self.tw_gridblocks.itemClicked.connect(
                self._grid_manager_item_click)
        self.tw_gridblocks.setContextMenuPolicy(
                QtCore.Qt.CustomContextMenu)
        self.tw_gridblocks.customContextMenuRequested.connect(
                self._grid_manager_context_menu)

    def _cont_manager_init(self):
        'initialization of contour manager dock'
        self.cont_view = contvis.ContManagerItemView(self.dw_contour_manager)
        self.dw_contour_manager.resize(200, 500)
        self.dw_contour_manager.setWidget(self.cont_view)

    def _bnd_manager_init(self):
        'initialization of boundary manager dock'
        self.bnd_view = btypes.BoundaryTypesView(self.dw_bnd_manager)
        self.dw_bnd_manager.resize(200, 500)
        self.dw_bnd_manager.setWidget(self.bnd_view)

    def _commands_history_init(self):
        'initialization of the commands history window'
        self.dw_history.setFloating(True)
        self.dw_history.resize(700, 300)
        self.dw_history.setShown(False)
        self.dw_history.setWidget(commandwin.HistoryDockFrame(
            globvars.Flows, self.dw_history))

    def _vtk_widget_init(self):
        ' initialization of the vtk renderer window'
        self.vtkWidget = drawarea.VTKWidget(self.f_vtk)
        layout = QHBoxLayout()
        layout.addWidget(self.vtkWidget)
        self.f_vtk.setLayout(layout)

    # --------- Flows Management
    def _save_flows(self):
        filename = QFileDialog.getSaveFileName(self,
            'Save to File', '', 'ComGrid projects (*.cgp);;All Files (*)')
        if (filename):
            globvars.Flows.xml_save(filename)

    def _open_flows(self):
        filename = QFileDialog.getOpenFileName(
            self, 'Open File', '', 'ComGrid projects (*.cgp);;All Files (*)')
        if (filename):
            globvars.Flows.xml_load(filename)

    def _undo_command(self):
        globvars.actual_flow().undo_prev()

    def _redo_command(self):
        globvars.actual_flow().exec_next()

    # --------- geometry actions
    def _add_rectangle_grid(self):
        dialog = dlgs.AddUnfRectGrid(self)
        if dialog.exec_():
            arg = dialog.ret_value()
            com = gridcom.AddUnfRectGrid(arg)
            globvars.actual_flow().exec_command(com)

    def _add_circular_grid(self):
        dialog = dlgs.AddUnfCircGrid(self)
        if dialog.exec_():
            r = dialog.ret_value()
            com = gridcom.AddUnfCircGrid(r)
            globvars.actual_flow().exec_command(com)

    def _add_ring_grid(self):
        dialog = dlgs.AddUnfRingGrid(self)
        if dialog.exec_():
            r = dialog.ret_value()
            com = gridcom.AddUnfRingGrid(r)
            globvars.actual_flow().exec_command(com)

    def _add_rectangle_cont(self):
        bnd_types = globvars.actual_data().get_bnd_types()
        dialog = dlgs.AddRectCont(bnd_types, self)
        if dialog.exec_():
            r = dialog.ret_value()
            com = contcom.AddRectCont(*r)
            globvars.actual_flow().exec_command(com)

    def _copy_grids(self):
        all_grids = globvars.actual_data().get_grid_names()
        used_grids = globvars.actual_data().get_checked_grid_names()
        all_conts = globvars.actual_data().get_contour_names()
        used_conts = globvars.actual_data().get_checked_contour_names()
        dialog = dlgs.CopyGrids(used_grids, all_grids,
                used_conts, all_conts, self)
        if dialog.exec_():
            srcnames, newnames, cnames, newcnames = dialog.ret_value()
            com = objcom.CopyGrid(srcnames, newnames, cnames, newcnames)
            globvars.actual_flow().exec_command(com)

    def _remove_grids(self):
        all_grids = globvars.actual_data().get_grid_names()
        used_grids = globvars.actual_data().get_checked_grid_names()
        all_conts = globvars.actual_data().get_contour_names()
        used_conts = globvars.actual_data().get_checked_contour_names()
        dialog = dlgs.RemoveGrids(used_grids, all_grids,
                used_conts, all_conts, self)
        if dialog.exec_():
            remnames, remcntnames = dialog.ret_value()
            com = objcom.RemoveGrid(remnames, remcntnames)
            globvars.actual_flow().exec_command(com)

    def _mov_rot(self):
        all_grids = globvars.actual_data().get_grid_names()
        used_grids = globvars.actual_data().get_checked_grid_names()
        all_conts = globvars.actual_data().get_contour_names()
        used_conts = globvars.actual_data().get_checked_contour_names()
        dialog = dlgs.MoveRotateGridsDlg(used_grids, all_grids,
                used_conts, all_conts, self)
        if dialog.exec_():
            names, cnames, mopt, ropt = dialog.ret_value()
            #1) rotate
            if abs(ropt[1]) > 1e-16:
                com = objcom.RotateGrids(ropt[0], ropt[1], names, cnames)
                globvars.actual_flow().exec_command(com)
            #2) move
            if abs(mopt[0]) > 1e-16 or abs(mopt[1]) > 1e-16:
                com = objcom.MoveGrids(mopt[0], mopt[1], names, cnames)
                globvars.actual_flow().exec_command(com)

    def _scale_grid(self):
        all_grids = globvars.actual_data().get_grid_names()
        used_grids = globvars.actual_data().get_checked_grid_names()
        all_conts = globvars.actual_data().get_contour_names()
        used_conts = globvars.actual_data().get_checked_contour_names()
        dialog = dlgs.ScaleGridsDlg(used_grids, all_grids,
                used_conts, all_conts, self)
        if dialog.exec_():
            names, cnames, rel_pnt, xpct, ypct = dialog.ret_value()
            if xpct != 100 or ypct != 100:
                com = objcom.ScaleGrids(rel_pnt, xpct, ypct, names, cnames)
                globvars.actual_flow().exec_command(com)

    # --------- tools
    def _unite_grids(self):
        all_grids = globvars.actual_data().get_grid_names()
        used_grids = globvars.actual_data().get_checked_grid_names()
        dialog = dlgs.UniteGrids(used_grids, all_grids, self)
        if dialog.exec_():
            name, pbnd, grd, bfs, den = dialog.ret_value()
            src = [gridcom.UniteOpts(n, b, d)
                    for n, b, d in zip(grd, bfs, den)]
            a = {'name': name, 'fix_bnd': pbnd}
            for i, s in enumerate(src):
                k = ''.join(['s', str(i)])
                a[k] = s
            com = gridcom.UniteGrids(**a)
            globvars.actual_flow().exec_command(com)

    def _set_bc(self):
        used_conts = globvars.actual_data().get_checked_any_contour_names()
        if len(used_conts) == 0:
            QtGui.QMessageBox.warning(None, "Warning", "Check contours first")
            return
        ret = globvars.actual_data().ask_for_new_contours_bnd(used_conts)
        if ret is not None:
            arg = []
            for c, e in zip(used_conts, ret):
                arg.append(contcom.BTypePicker(c, e))
            com = contcom.SetBTypeToContour(arg)
            globvars.actual_flow().exec_command(com)

    # --------- Grid Manager
    def clear_grid_manager(self):
        ' clears QTreeWidgetItem '
        self.tw_gridblocks.clear()
        self.tw_gridblocks.setColumnCount(5)

    def add_grid_manager_item(self, item):
        ' adds QTreeWidgetItem to the Grid Manager '
        #add item
        self.tw_gridblocks.invisibleRootItem().addChild(item)
        #resize column width
        for i in range(self.tw_gridblocks.columnCount()):
            self.tw_gridblocks.resizeColumnToContents(i)

    def _grid_manager_item_click(self, item, col):
        item.clicked(col)

    def _grid_manager_context_menu(self, pnt):
        col = self.tw_gridblocks.indexAt(pnt).column()
        if col == 0:
            item = self.tw_gridblocks.itemAt(pnt)
            item.show_context(pnt)

    # ---------- Interfaces
    def actual_flow_changed(self, flow_name):
        "subscription from command.FlowCollection"
        self.setWindowTitle("HybMesh: %s" % flow_name)
