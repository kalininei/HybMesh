#!/usr/bin/env python

from PyQt4.QtGui import (QMainWindow,
    QFileDialog, QHBoxLayout)
import qtui.ui_ComGridMain

import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

import gridcom
import dlgs
import globvars


class VTKWidget(QVTKRenderWindowInteractor):
    def __init__(self, parent=None):
        super(VTKWidget, self).__init__(parent)
        #renderer
        self.ren = vtk.vtkRenderer()
        self.GetRenderWindow().AddRenderer(self.ren)
        self.ren.SetBackground(0, 0.3, 0.3)
        #2d interactor style
        self.iren().SetInteractorStyle(vtk.vtkInteractorStyleImage())

    def iren(self):
        return self.GetRenderWindow().GetInteractor()

    def _get_actors(self):
        ret = []
        it = self.ren.GetActors().NewIterator()
        while not it.IsDoneWithTraversal():
            ret.append(it.GetCurrentObject())
            it.GoToNextItem()
        return ret

    def set_actor_list(self, actlist):
        ' sets actors for visualization '
        #removes all actors which are not in actlist
        for a in self._get_actors():
            if a not in actlist:
                self.ren.RemoveActor(a)
        #add actlist
        #no repeats due to AddActor implementation
        map(self.ren.AddActor, actlist)


class MainWindow(QMainWindow, qtui.ui_ComGridMain.Ui_MainWindow):
    ' Main program window '
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)
        self._connect_actions()
        self._grid_manager_init()
        self._vtk_widget_init()

    # --------- constructor procedures
    def _connect_actions(self):
        ' connects actions for main window'
        #file
        self.act_open.triggered.connect(self._open_flows)
        self.act_save.triggered.connect(self._save_flows)
        #edit
        self.act_undo.triggered.connect(self._undo_command)
        self.act_redo.triggered.connect(self._redo_command)
        #view
        self.act_comhist.triggered.connect(self._show_history)
        #geometry
        self.act_unf_rect.triggered.connect(self._add_rectangle_grid)
        self.act_unf_circ.triggered.connect(self._add_circular_grid)
        self.act_unf_ring.triggered.connect(self._add_ring_grid)
        #tools
        self.act_unite_grids.triggered.connect(self._unite_grids)
        #transform
        self.act_movrot.triggered.connect(self._mov_rot)
        self.act_scale.triggered.connect(self._scale_grid)

    def _grid_manager_init(self):
        ' initialization of the grid manager '
        self.tw_gridblocks.itemClicked.connect(
                self._grid_manager_item_click)

    def _vtk_widget_init(self):
        ' initialization of the vtk renderer window'
        self.vtkWidget = VTKWidget(self.f_vtk)
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

    def _show_history(self):
        globvars.Flows.show_manager()

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

    # ---------- transform
    def _mov_rot(self):
        all_grids = globvars.actual_data().get_grid_names()
        used_grids = globvars.actual_data().get_checked_grid_names()
        dialog = dlgs.MoveRotateGridsDlg(used_grids, all_grids, self)
        if dialog.exec_():
            names, mopt, ropt = dialog.ret_value()
            #1) rotate
            if abs(ropt[1]) > 1e-16:
                com = gridcom.RotateGrids(ropt[0], ropt[1], names)
                globvars.actual_flow().exec_command(com)
            #2) move
            if abs(mopt[0]) > 1e-16 or abs(mopt[1]) > 1e-16:
                com = gridcom.MoveGrids(mopt[0], mopt[1], names)
                globvars.actual_flow().exec_command(com)

    def _scale_grid(self):
        all_grids = globvars.actual_data().get_grid_names()
        used_grids = globvars.actual_data().get_checked_grid_names()
        dialog = dlgs.ScaleGridsDlg(used_grids, all_grids, self)
        if dialog.exec_():
            names, rel_pnt, xpct, ypct = dialog.ret_value()
            if xpct != 100 or ypct != 100:
                com = gridcom.ScaleGrids(rel_pnt, xpct, ypct, names)
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

    # ---------- Grid Viewer
    def set_grid_actors(self, actors):
        self.vtkWidget.set_actor_list(actors)
