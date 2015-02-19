#!/usr/bin/env python

from PyQt4 import QtCore
from PyQt4.QtGui import (QDialog, QDialogButtonBox, QGridLayout,
    QLabel, QLineEdit, QMessageBox, QProgressBar, QVBoxLayout)
import bgeom
import qtui.ui_UniteGridsDlg
import qtui.ui_AddUnfCircDlg
import qtui.ui_MoveRotateDlg
import qtui.ui_GridViewOpt


class _BackGroundWorkerCB(QtCore.QThread):
    """ Background procedure with callback caller
        for call from ProgressProcedureDlg
    """
    def __init__(self, emitter, func, args):
        """ emitter - QObject that emits signal,
            func(*args, callback_fun) - target procedure,
            callback function should be declared as
                int callback_fun(QString BaseName, QString SubName,
                    double proc1, double proc2)
            it should return 1 for cancellation requiry and 0 otherwise

        """
        super(_BackGroundWorkerCB, self).__init__()
        self.emitter = emitter
        self.func = func
        self.args = args + (self._get_callback(),)
        self.proceed, self._result = True, None

    def run(self):
        self.proceed = True
        self._result = self.func(*self.args)

    def _emit(self, n1, n2, p1, p2):
        self.emitter.emit(QtCore.SIGNAL(
            "fill_cb_form(QString, QString, double, double)"), n1, n2, p1, p2)

    def _get_callback(self):
        import ctypes as ct

        def cb(n1, n2, p1, p2):
            self._emit(n1, n2, p1, p2)
            return 0 if self.proceed else 1

        #if target function is a c function then convert callback to
        #a c function pointer
        if isinstance(self.func, ct._CFuncPtr):
            cbfunc = ct.CFUNCTYPE(ct.c_int, ct.c_char_p, ct.c_char_p,
                   ct.c_double, ct.c_double)
            cb2 = cbfunc(cb)
            return cb2
        else:
            return cb


class ProgressProcedureDlg(QDialog):
    """ ProgressBar/Cancel dialog window which wraps time consuming
        procedure calls.
    """

    def __init__(self, func, args, parent=None):
        """ func(*args, callback_fun) - target procedure,
            callback_fun should be declared as
                int callback_fun(QString BaseName, QString SubName,
                    double proc1, double proc2)
            it should return 1 for cancellation requiry and 0 otherwise
        """
        flags = QtCore.Qt.Dialog \
                | QtCore.Qt.CustomizeWindowHint \
                | QtCore.Qt.WindowTitleHint
        super(ProgressProcedureDlg, self).__init__(parent, flags)

        #design
        self._label1 = QLabel()
        self._label2 = QLabel()
        self._progbar1 = QProgressBar()
        self._progbar2 = QProgressBar()
        self._buttonbox = QDialogButtonBox(QDialogButtonBox.Cancel)
        self._buttonbox.setCenterButtons(True)
        self._buttonbox.rejected.connect(self._cancel_pressed)
        layout = QVBoxLayout()
        layout.addWidget(self._label1)
        layout.addWidget(self._progbar1)
        layout.addWidget(self._label2)
        layout.addWidget(self._progbar2)
        layout.addWidget(self._buttonbox)
        self.setLayout(layout)
        self.setFixedSize(400, 150)
        self.setModal(True)

        #worker
        emitter = QtCore.QObject()
        self.connect(emitter, QtCore.SIGNAL(
            "fill_cb_form(QString, QString, double, double)"), self._fill)

        self._worker = _BackGroundWorkerCB(emitter, func, args)
        self._worker.finished.connect(self._fin)
        self._worker.terminated.connect(self._fin)

    def exec_(self):
        self.show()  # show first to avoid delay
        self._worker.start()
        return super(ProgressProcedureDlg, self).exec_()

    def _fill(self, n1, n2, p1, p2):
        self.setWindowTitle(n1)
        self._label1.setText(n1)
        self._label2.setText(n2)
        self._progbar1.setValue(100 * p1)
        self._progbar2.setValue(100 * p2)

    def _cancel_pressed(self):
        self._worker.proceed = False

    def _fin(self):
        self._result = self._worker._result
        self.close()

    def get_result(self):
        return self._result


class AddUnfRectGrid(QDialog):
    ' Add Uniform Rectangular Grid dialog '
    def __init__(self, parent=None):
        super(AddUnfRectGrid, self).__init__(parent)
        #buttons
        buttonbox = QDialogButtonBox(
                QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttonbox.accepted.connect(self.accept)
        buttonbox.rejected.connect(self.reject)
        #edits
        self.p0xedit = QLineEdit()
        self.p0xedit.setText("0")
        self.p0yedit = QLineEdit()
        self.p0yedit.setText("0")
        self.p1xedit = QLineEdit()
        self.p1xedit.setText("1")
        self.p1yedit = QLineEdit()
        self.p1yedit.setText("1")
        self.nmedit = QLineEdit()
        self.nmedit.setText("")
        self.nxedit = QLineEdit()
        self.nxedit.setText("10")
        self.nyedit = QLineEdit()
        self.nyedit.setText("10")
        #Layout
        layout = QGridLayout()

        #Grid name
        layout.addWidget(QLabel("Grid name"), 3, 0, 1, 4)
        layout.addWidget(self.nmedit, 4, 0, 1, 4)

        #Grid partition
        layout.addWidget(QLabel("Grid partition"), 7, 0, 1, 4)
        layout.addWidget(QLabel("Nx"), 8, 0)
        layout.addWidget(QLabel("Ny"), 8, 2)
        layout.addWidget(self.nxedit, 8, 1)
        layout.addWidget(self.nyedit, 8, 3)

        #left bottom
        layout.addWidget(QLabel("bottom left:"), 10, 0, 1, 4)
        layout.addWidget(QLabel("X"), 20, 0)
        layout.addWidget(self.p0xedit, 20, 1)
        layout.addWidget(QLabel("Y"), 20, 2)
        layout.addWidget(self.p0yedit, 20, 3)
        #right top
        layout.addWidget(QLabel("top right"), 30, 0, 1, 4)
        layout.addWidget(QLabel("X"), 40, 0)
        layout.addWidget(self.p1xedit, 40, 1)
        layout.addWidget(QLabel("Y"), 40, 2)
        layout.addWidget(self.p1yedit, 40, 3)

        layout.addWidget(buttonbox, 100, 0, 1, 4)
        self.setLayout(layout)

    def accept(self):
        try:
            self.ret_value()
            super(AddUnfRectGrid, self).accept()
        except Exception:
            QMessageBox.warning(self, "Warning", "Invalid Input")

    def ret_value(self):
        ' -> (p0, p1, Nx, Ny, GridName) '
        p0 = bgeom.Point2(float(self.p0xedit.text()),
                float(self.p0yedit.text()))
        p1 = bgeom.Point2(float(self.p1xedit.text()),
                float(self.p1yedit.text()))
        Nx, Ny = int(self.nxedit.text()), int(self.nyedit.text())
        name = str(self.nmedit.text())
        #input data check
        if Nx <= 0 or Ny <= 0 or p0.x >= p1.x or p0.y >= p1.y:
            raise Exception
        return (p0, p1, Nx, Ny, name)


class AddUnfCircGrid(QDialog, qtui.ui_AddUnfCircDlg.Ui_add_unf_circ):
    ' Add Uniform circular grid dialog '

    def __init__(self, parent=None):
        super(AddUnfCircGrid, self).__init__(parent)
        self.setupUi(self)

    def ret_value(self):
        '-> (pc, rad, Na, Nr, ref_coef, trian_center, GridName) '
        pc = bgeom.Point2(float(self.ed_x.text()),
                float(self.ed_y.text()))
        rad = float(self.ed_rad.text())
        Na, Nr = int(self.ed_na.text()), int(self.ed_nr.text())
        ref_coef = float(self.ed_coef.text())
        trian_center = self.cb_trian.isChecked()
        name = str(self.ed_name.text())
        return pc, rad, Na, Nr, ref_coef, trian_center, name

    def accept(self):
        try:
            self.ret_value()
            super(AddUnfCircGrid, self).accept()
        except Exception:
            QMessageBox.warning(self, "Warning", "Invalid Input")


class UniteGrids(QDialog, qtui.ui_UniteGridsDlg.Ui_Dialog):
    ' unite grids dialog window '

    def __init__(self, act_grids, inact_grids):
        super(UniteGrids, self).__init__()
        self._act = act_grids
        self._inact = inact_grids
        self.setupUi(self)
        self._set_lists()
        self.bt_left.clicked.connect(self._bleft)
        self.bt_right.clicked.connect(self._bright)

    def _set_lists(self):
        self.lst_used.clear()
        self.lst_unused.clear()
        map(self.lst_used.addItem, self._act)
        map(self.lst_unused.addItem, self._inact)

    def _bleft(self):
        for it in self.lst_unused.selectedItems():
            self._act.append(it.text())
            self._inact.remove(it.text())
        self._set_lists()

    def _bright(self):
        for it in self.lst_used.selectedItems():
            self._inact.append(it.text())
            self._act.remove(it.text())
        self._set_lists()

    def accept(self):
        ' overridden accept with input check '
        try:
            self.ret_value()
            super(UniteGrids, self).accept()
        except Exception:
            QMessageBox.warning(self, "Warning", "Invalid Input")

    def ret_value(self):
        ' -> (GridName, [source grids], [buffer sizes], [densities]) '
        nm = str(self.ed_name.text())
        src, buf, den = [], [], []
        for i in range(self.lst_used.count()):
            it = self.lst_used.item(i)
            src.append(str(it.text()))
            buf.append(float(self.ed_buf.text()))
            den.append(self.sld_dens.value())
        if len(src) < 2:
            #not enough grids to unite
            raise Exception()
        return nm, src, buf, den


class MoveRotateGridsDlg(QDialog, qtui.ui_MoveRotateDlg.Ui_Dialog):
    ' Move Rotate grids group'

    def __init__(self, act_grids, inact_grids, parent=None):
        super(MoveRotateGridsDlg, self).__init__(parent)
        self.setupUi(self)
        self.transfer_list.set_lists(act_grids, inact_grids)

    def accept(self):
        ' overridden accept with input check '
        try:
            self.ret_value()
            super(MoveRotateGridsDlg, self).accept()
        except Exception:
            QMessageBox.warning(self, "Warning", "Invalid Input")

    def ret_value(self):
        ' [source grids], [move_x, move_y], [rot_x0, rot_y0, rot_angle] '
        names = self.transfer_list.get_left_list()
        mx, my = float(self.ed_movx.text()), float(self.ed_movy.text())
        rx0, ry0 = float(self.ed_rotx0.text()), float(self.ed_roty0.text())
        an = float(self.ed_rotan.text())
        if len(names) == 0:
            raise Exception()
        return names, [mx, my], [rx0, ry0, an]


class GridViewOpt(QDialog, qtui.ui_GridViewOpt.Ui_Dialog):
    ' Grid view options '

    def __init__(self, opts, parent=None):
        super(GridViewOpt, self).__init__(parent)
        self.setupUi(self)
        self.apply_opts(opts)

    def apply_opts(self, opts):
        #TODO
        pass


if __name__ == "__main__":
    import sys
    from PyQt4.QtGui import QApplication
    app = QApplication(sys.argv)
    d = AddUnfRectGrid()
    d.exec_()

    s = ";".join(map(str, d.retValue()))
    print s
