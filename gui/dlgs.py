#!/usr/bin/env python

import copy
from PyQt4 import QtCore, QtGui
import bgeom
import globvars
import optview
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


class ProgressProcedureDlg(QtGui.QDialog):
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
        self._label1 = QtGui.QLabel()
        self._label2 = QtGui.QLabel()
        self._progbar1 = QtGui.QProgressBar()
        self._progbar2 = QtGui.QProgressBar()
        self._buttonbox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel)
        self._buttonbox.setCenterButtons(True)
        self._buttonbox.rejected.connect(self._cancel_pressed)
        layout = QtGui.QVBoxLayout()
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


class SimpleAbstractDialog(QtGui.QDialog):
    "Abstract dialog for option set"
    class _OData(object):
        pass

    def odata(self):
        """returns options struct child class singleton which stores
        last dialog execution"""
        if not hasattr(self, "_odata"):
            setattr(self.__class__, "_odata", SimpleAbstractDialog._OData())
            self._default_odata(self._odata)
        return self._odata

    def __init__(self, parent=None):
        super(SimpleAbstractDialog, self).__init__(parent)
        oview = optview.OptionsView(self.olist())
        oview.is_active_delegate(self._active_entries)
        buttonbox = QtGui.QDialogButtonBox(
                QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel)
        buttonbox.accepted.connect(self.accept)
        buttonbox.rejected.connect(self.reject)
        layout = QtGui.QVBoxLayout(self)
        layout.addWidget(oview)
        layout.addWidget(buttonbox)

    def accept(self):
        "check errors and invoke parent accept"
        try:
            self.check_input()
            super(SimpleAbstractDialog, self).accept()
        except Exception as e:
            QtGui.QMessageBox.warning(self, "Warning",
                    "Invalid Input: %s" % str(e))

    #functions for overriding
    def _default_odata(self, obj):
        "fills options struct obj with default values"
        raise NotImplementedError

    def olist(self):
        "-> optview.OptionsList"
        raise NotImplementedError

    def ret_value(self):
        "-> return value from option struct"
        raise NotImplementedError

    def check_input(self):
        "throws Exception if self.odata() has invalid fields"
        pass

    def _active_entries(self, entry):
        """ (optview.OptionEntry entry) -> bool
            return False for non-active entries
        """
        return True


class AddUnfRectGrid(SimpleAbstractDialog):
    'Add Uniform Rectangular Grid dialog '

    def __init__(self, parent=None):
        super(AddUnfRectGrid, self).__init__(parent)
        self.resize(300, 300)
        self.setWindowTitle("Build uniform rectangular grid")

    def _default_odata(self, obj):
        "-> options struct with default values"
        obj.name = "RectGrid1"
        obj.p0, obj.p1 = bgeom.Point2(0.0, 0.0), bgeom.Point2(1.0, 1.0)
        obj.nx, obj.ny = 10, 10

    def olist(self):
        "-> optview.OptionsList"
        return optview.OptionsList([("Basic", "Grid name",
                optview.SimpleOptionEntry(self.odata(), "name")),
            ("Geometry", "Bottom left",
                optview.XYOptionEntry(self.odata(), "p0")),
            ("Geometry", "Top right",
                optview.XYOptionEntry(self.odata(), "p1")),
            ("Partition", "Nx",
                optview.BoundedIntOptionEntry(self.odata(), "nx", 1, 1e6)),
            ("Partition", "Ny",
                optview.BoundedIntOptionEntry(self.odata(), "ny", 1, 1e6))])

    def ret_value(self):
        ' -> {p0, p1, nx, ny, name} '
        od = copy.deepcopy(self.odata())
        return {"p0": od.p0, "p1": od.p1, "nx": od.nx, "ny": od.ny,
                "name": od.name}

    def check_input(self):
        if self.odata().p0.x == self.odata().p1.x or \
                self.odata().p0.y == self.odata().p1.y:
            raise Exception("Zero area polygon")
        if self.odata().p0.x >= self.odata().p1.x or \
                self.odata().p0.y >= self.odata().p1.y:
            raise Exception("Invalid points order")


class AddUnfCircGrid(SimpleAbstractDialog):
    ' Add Uniform circular grid dialog '

    def __init__(self, parent=None):
        super(AddUnfCircGrid, self).__init__(parent)
        self.resize(400, 400)
        self.setWindowTitle("Build uniform circular grid")

    def _default_odata(self, obj):
        "-> options struct with default values"
        obj.name = "CircGrid1"
        obj.pc, obj.rad = bgeom.Point2(0.0, 0.0), 3.0
        obj.na, obj.nr = 10, 4
        obj.ctrian, obj.coef = True, 1.2

    def olist(self):
        "-> optview.OptionsList"
        return optview.OptionsList([("Basic", "Grid name",
                optview.SimpleOptionEntry(self.odata(), "name")),
            ("Geometry", "Center point",
                optview.XYOptionEntry(self.odata(), "pc")),
            ("Geometry", "Radius",
                optview.SimpleOptionEntry(self.odata(), "rad")),
            ("Partition", "Radius partition",
                optview.BoundedIntOptionEntry(self.odata(), "nr", 1, 1e6)),
            ("Partition", "Arch partition",
                optview.BoundedIntOptionEntry(self.odata(), "na", 3, 1e6)),
            ("Partition", "Refinement coef.",
                optview.SimpleOptionEntry(self.odata(), "coef")),
            ("Partition", "Triangulate center cell",
                optview.BoolOptionEntry(self.odata(), "ctrian"))])

    def check_input(self):
        if self.odata().rad <= 0:
            raise Exception("Invalid radius data")

    def ret_value(self):
        '-> {pc, rad, Na, Nr, coef, is_trian, name}'
        od = copy.deepcopy(self.odata())
        return {'p0': od.pc, 'rad': od.rad, 'na': od.na, 'nr': od.nr,
                'coef': od.coef, 'is_trian': od.ctrian, 'name': od.name}


class AddUnfRingGrid(SimpleAbstractDialog):
    ' Add Uniform circular ring dialog '

    def __init__(self, parent=None):
        super(AddUnfRingGrid, self).__init__(parent)
        self.resize(400, 400)
        self.setWindowTitle("Build uniform ring grid")

    def _default_odata(self, obj):
        "-> options struct with default values"
        obj.name = "RingGrid1"
        obj.pc, obj.irad, obj.orad = bgeom.Point2(0.0, 0.0), 1.0, 4.0
        obj.na, obj.nr = 10, 4
        obj.coef = 1.2

    def olist(self):
        "-> optview.OptionsList"
        return optview.OptionsList([("Basic", "Grid name",
                optview.SimpleOptionEntry(self.odata(), "name")),
            ("Geometry", "Center point",
                optview.XYOptionEntry(self.odata(), "pc")),
            ("Geometry", "Inner Radius",
                optview.SimpleOptionEntry(self.odata(), "irad")),
            ("Geometry", "Outer Radius",
                optview.SimpleOptionEntry(self.odata(), "orad")),
            ("Partition", "Radius partition",
                optview.BoundedIntOptionEntry(self.odata(), "nr", 1, 1e6)),
            ("Partition", "Arch partition",
                optview.BoundedIntOptionEntry(self.odata(), "na", 3, 1e6)),
            ("Partition", "Refinement coefficient",
                optview.SimpleOptionEntry(self.odata(), "coef"))])

    def ret_value(self):
        '-> {pc, radinner, radouter, na, nr, coef, name}'
        od = copy.deepcopy(self.odata())
        return {'p0': od.pc, 'radinner': od.irad, 'radouter': od.orad,
                'na': od.na, 'nr': od.nr, 'coef': od.coef, 'name': od.name}

    def check_input(self):
        if self.odata().irad <= 0 or self.odata().orad <= 0:
            raise Exception("invalid radius data")
        if self.odata().irad >= self.odata().orad:
            raise Exception("inner radius is greater then outer")


class UniteGrids(SimpleAbstractDialog):
    ' unite grids dialog window '

    def __init__(self, used_grids, all_grids, parent=None):
        self.all_grids = all_grids
        self.odata().grds = used_grids[:]
        super(UniteGrids, self).__init__(parent)
        self.resize(400, 300)
        self.setWindowTitle("Unite grids")

    def _default_odata(self, obj):
        "-> options struct with default values"
        obj.name = "UnitedGrid1"
        obj.pbnd = False
        obj.buff = 0.3
        obj.den = 7
        obj.grds = []

    def olist(self):
        "-> optview.OptionsList"
        return optview.OptionsList([("Basic", "Grid name",
                optview.SimpleOptionEntry(self.odata(), "name")),
            ("Construction", "Preserve boundary nodes",
                optview.BoolOptionEntry(self.odata(), "pbnd")),
            ("Construction", "Buffer size",
                optview.SimpleOptionEntry(self.odata(), "buff")),
            ("Construction", "Density",
                optview.BoundedIntOptionEntry(self.odata(), "den", 0, 10)),
            ("Grids", "Grids list", optview.MultipleChoiceOptionEntry(
                self.odata(), "grds", self.all_grids))])

    def check_input(self):
        if self.odata().buff < 0:
            raise Exception("Invalid buffer size")
        if len(self.odata().grds) < 2:
            raise Exception("not enough grids")

    def ret_value(self):
        """ -> (GridName, preserve_bnd,
            [source grids], [buffer sizes], [densities])
        """
        od = copy.deepcopy(self.odata())
        buff = [od.buff] * len(od.grds)
        den = [od.den] * len(od.grds)
        return od.name, od.pbnd, od.grds, buff, den


class MoveRotateGridsDlg(SimpleAbstractDialog):
    ' Move Rotate grids group'

    def __init__(self, used_grids, all_grids,
            used_conts, all_conts, parent=None):
        self.odata().grds = used_grids
        self.odata().cnts = used_conts
        self.all_grids = all_grids
        self.all_conts = all_conts
        super(MoveRotateGridsDlg, self).__init__(parent)
        self.resize(400, 300)
        self.setWindowTitle("Move/Rotate objects")

    def _default_odata(self, obj):
        "-> options struct with default values"
        obj.do_move, obj.do_rotate = True, True
        obj.dx, obj.dy = 0.0, 0.0
        obj.rotp = bgeom.Point2(0.0, 0.0)
        obj.rotangle = 0.0
        obj.grds = []
        obj.cnts = []

    def olist(self):
        "-> optview.OptionsList"
        return optview.OptionsList([("Move", "Move grids",
                optview.BoolOptionEntry(self.odata(), "do_move")),
            ("Move", "X shift",
                optview.SimpleOptionEntry(self.odata(), "dx")),
            ("Move", "Y shift",
                optview.SimpleOptionEntry(self.odata(), "dy")),
            ("Rotate", "Rotate grids",
                optview.BoolOptionEntry(self.odata(), "do_rotate")),
            ("Rotate", "Reference Point",
                optview.XYOptionEntry(self.odata(), "rotp")),
            ("Rotate", "Angle (deg)",
                optview.SimpleOptionEntry(self.odata(), "rotangle")),
            ("Objects", "Grids list", optview.MultipleChoiceOptionEntry(
                self.odata(), "grds", self.all_grids)),
            ("Objects", "Contours list", optview.MultipleChoiceOptionEntry(
                self.odata(), "cnts", self.all_conts))])

    def _active_entries(self, entry):
        "return False for non-active entries"
        if entry.member_name in ["dx", "dy"]:
            return entry.data.do_move
        elif entry.member_name in ["rotp", "rotangle"]:
            return entry.data.do_rotate
        elif entry.member_name in ["grds"]:
            return entry.data.do_move or entry.data.do_rotate
        else:
            return True

    def check_input(self):
        if (len(self.odata().grds) < 1 and len(self.odata().cnts) < 1) or \
                (not self.odata().do_move and not self.odata().do_rotate):
            raise Exception("nothing to do")

    def ret_value(self):
        """->[source grids], [source conts],
             [move_x, move_y], [rot_point, rot_angle]
        """
        od = copy.deepcopy(self.odata())
        if not od.do_move:
            od.dx = od.dy = 0
        if not od.do_rotate:
            od.rotangle = 0
        return od.grds, od.cnts, [od.dx, od.dy], [od.rotp, od.rotangle]


class ScaleGridsDlg(SimpleAbstractDialog):
    'Scale objects group'

    def __init__(self, used_grids, all_grids,
            used_cnts, all_cnts, parent=None):
        self.odata().grds = used_grids
        self.odata().cnts = used_cnts
        self.all_grids = all_grids
        self.all_cnts = all_cnts
        super(ScaleGridsDlg, self).__init__(parent)
        self.resize(400, 300)
        self.setWindowTitle("Scale objects")

    def _default_odata(self, obj):
        "-> options struct with default values"
        obj.width = 100.0
        obj.height = 100.0
        obj.units = "%"
        obj.preserve_ratio = True
        obj.point = "center"
        obj.grds = []
        obj.cnts = []

    def olist(self):
        "-> optview.OptionsList"
        a_pnt = ["bottom left", "center", "top right"]
        a_units = ["%", "meter"]
        return optview.OptionsList([("Scale", "Preserve ratio",
                optview.BoolOptionEntry(self.odata(), "preserve_ratio")),
            ("Scale", "Units", optview.SingleChoiceOptionEntry(
                self.odata(), "units", a_units)),
            ("Scale", "Width",
                optview.SimpleOptionEntry(self.odata(), "width")),
            ("Scale", "Height",
                optview.SimpleOptionEntry(self.odata(), "height")),
            ("Reference", "Reference Point", optview.SingleChoiceOptionEntry(
                self.odata(), "point", a_pnt)),
            ("Objects", "Grids list", optview.MultipleChoiceOptionEntry(
                self.odata(), "grds", self.all_grids)),
            ("Objects", "Contours list", optview.MultipleChoiceOptionEntry(
                self.odata(), "cnts", self.all_cnts))])

    def _active_entries(self, entry):
        "return False for non-active entries"
        if entry.member_name in ["height"]:
            return not entry.data.preserve_ratio
        else:
            return True

    def check_input(self):
        if (len(self.odata().grds) < 1 and len(self.odata().cnts) < 1):
            raise Exception("nothing to do")
        if self.odata().width <= 0 or \
                (not self.odata().preserve_ratio and self.odata().height <= 0):
            raise Exception("invalid scales")

    def ret_value(self):
        "-> names, contnames, rel_pnt, xpct, ypct"
        od = copy.deepcopy(self.odata())

        #builing bounding box
        x, y = [], []
        for n in od.grds:
            p0, p1 = globvars.actual_data().grids2[n].bounding_box()
            x += [p0.x, p1.x]
            y += [p0.y, p1.y]
        for n in od.cnts:
            p0, p1 = globvars.actual_data().contours2[n].bounding_box()
            x += [p0.x, p1.x]
            y += [p0.y, p1.y]
        p0, p1 = bgeom.Point2(min(x), min(y)), bgeom.Point2(max(x), max(y))

        #relative point
        if od.point == "bottom left":
            rel_pnt = p0
        elif od.point == "top right":
            rel_pnt = p1
        else:
            rel_pnt = bgeom.Point2((p0.x + p1.x) / 2.0, (p0.y + p1.y) / 2.0)

        #scaling units to %
        sx, sy = od.width, od.height
        if od.units == "meter":
            sx = sx / (p1.x - p0.x) * 100
            if (not od.preserve_ratio):
                sy = sy / (p1.y - p0.y) * 100
        if (od.preserve_ratio):
            sy = sx

        return od.grds, od.cnts, rel_pnt, sx, sy


class CopyGrids(SimpleAbstractDialog):
    def __init__(self, used_grids, all_grids,
            used_cnts, all_cnts, parent=None):
        self.odata().grds = used_grids
        self.odata().cnts = used_cnts
        self.all_grids = all_grids
        self.all_cnts = all_cnts
        super(CopyGrids, self).__init__(parent)
        self.resize(400, 300)
        self.setWindowTitle("Copy objects")

    def _default_odata(self, obj):
        "fills options struct obj with default values"
        obj.prefix = "Copy"
        obj.postfix = ""
        obj.grds = []
        obj.cnts = []

    def olist(self):
        "-> optview.OptionsList"
        return optview.OptionsList([("Names", "Grid name prefix",
                optview.SimpleOptionEntry(self.odata(), "prefix")),
            ("Names", "Grid name postfix",
                optview.SimpleOptionEntry(self.odata(), "postfix")),
            ("Objects", "Grids list", optview.MultipleChoiceOptionEntry(
                self.odata(), "grds", self.all_grids)),
            ("Objects", "Contours list", optview.MultipleChoiceOptionEntry(
                self.odata(), "cnts", self.all_cnts))])

    def ret_value(self):
        "-> (srcnames, newnames, srccontnames, newcontnames)"
        names = self.odata().grds[:]
        contnames = self.odata().cnts[:]
        newnames, cnewnames = [], []
        for n1 in names:
            newnames.append(self.odata().prefix + n1 + self.odata().postfix)
        for n1 in contnames:
            cnewnames.append(self.odata().prefix + n1 + self.odata().postfix)
        return names, newnames, contnames, cnewnames

    def check_input(self):
        "throws Exception if self.odata() has invalid fields"
        if len(self.odata().grds) < 1 and len(self.odata().cnts) < 1:
            raise Exception("No source objects")


class RemoveGrids(SimpleAbstractDialog):
    def __init__(self, used_grids, all_grids,
            used_conts, all_conts, parent=None):
        self.odata().grds = used_grids
        self.odata().cnts = used_conts
        self.all_grids = all_grids
        self.all_conts = all_conts
        super(RemoveGrids, self).__init__(parent)
        self.resize(400, 300)
        self.setWindowTitle("Remove objects")

    def _default_odata(self, obj):
        "fills options struct obj with default values"
        obj.grds = []
        obj.cnts = []

    def olist(self):
        "-> optview.OptionsList"
        from optview import MultipleChoiceOptionEntry as Mce
        return optview.OptionsList([
            ("Remove Objects", "Grids list", Mce(
                self.odata(), "grds", self.all_grids)),
            ("Remove Objects", "Contours list", Mce(
                self.odata(), "cnts", self.all_conts)),
        ])

    def ret_value(self):
        "-> [names], [contour names]"
        return self.odata().grds[:], self.odata().cnts[:]

    def check_input(self):
        "throws Exception if self.odata() has invalid fields"
        if len(self.odata().grds) < 1 and \
                len(self.odata().cnts) < 1:
            raise Exception("No objects to remove")


class AddRectCont(SimpleAbstractDialog):
    'Add Uniform Rectangular Grid dialog '

    def __init__(self, bnd_types, parent=None):
        self.bnd_types = bnd_types
        super(AddRectCont, self).__init__(parent)
        self.resize(300, 300)
        self.setWindowTitle("Build rectangular contour")

    def _default_odata(self, obj):
        "-> options struct with default values"
        obj.name = "RectCont1"
        obj.p0, obj.p1 = bgeom.Point2(0.0, 0.0), bgeom.Point2(1.0, 1.0)
        obj.bt_separate = False
        obj.bt_b = obj.bt_r = obj.bt_t = obj.bt_l = "0: None"

    def olist(self):
        "-> optview.OptionsList"
        self._default_odata(self.odata())
        names = self.bnd_types.get_names()
        caps = []
        for n in names:
            b = self.bnd_types.get(name=n)
            caps.append("%i: %s" % (b.index, b.name))
        return optview.OptionsList([("Basic", "Contour name",
                optview.SimpleOptionEntry(self.odata(), "name")),
            ("Geometry", "Bottom left",
                optview.XYOptionEntry(self.odata(), "p0")),
            ("Geometry", "Top right",
                optview.XYOptionEntry(self.odata(), "p1")),
            ("Boundary types", "Separate types",
                optview.BoolOptionEntry(self.odata(), "bt_separate")),
            ("Boundary types", "Bottom",
                optview.SingleChoiceOptionEntry(self.odata(), "bt_b", caps)),
            ("Boundary types", "Right",
                optview.SingleChoiceOptionEntry(self.odata(), "bt_r", caps)),
            ("Boundary types", "Top",
                optview.SingleChoiceOptionEntry(self.odata(), "bt_t", caps)),
            ("Boundary types", "Left",
                optview.SingleChoiceOptionEntry(self.odata(), "bt_l", caps)),
        ])

    def ret_value(self):
        ' -> (name, p0, p1, [bot, right, top, left bnd int indicies])'
        def to_b(s):
            r = self.odata().__dict__[s]
            r = r.split(None, 2)[0][:-1]
            return int(r)

        if self.odata().bt_separate:
            tp = map(to_b, ["bt_b", "bt_r", "bt_t", "bt_l"])
        else:
            tp = [to_b("bt_b")] * 4
        return self.odata().name, self.odata().p0, self.odata().p1, tp

    def check_input(self):
        if self.odata().p0.x == self.odata().p1.x or \
                self.odata().p0.y == self.odata().p1.y:
            raise Exception("Zero area polygon")
        if self.odata().p0.x >= self.odata().p1.x or \
                self.odata().p0.y >= self.odata().p1.y:
            raise Exception("Invalid points order")

    def _active_entries(self, entry):
        """ (optview.OptionEntry entry) -> bool
            return False for non-active entries
        """
        if entry.member_name in ["bt_r", "bt_t", "bt_l"]:
            return entry.data.bt_separate
        else:
            return True


class EditBoundaryType(SimpleAbstractDialog):
    'Add Uniform Rectangular Grid dialog '

    def __init__(self, bnd_types, bt, parent=None):
        self._fill_odata(bnd_types, bt, self.odata())
        #if zero boundary => only color dialog is active
        if bt is not None and bnd_types.get(name=bt.name).index == 0:
            self._ind0 = 0
        else:
            self._ind0 = 1
        super(EditBoundaryType, self).__init__(parent)
        self.resize(300, 300)
        self.setWindowTitle("Edit boundary type")

    def _default_odata(self, obj):
        "-> options struct with default values"
        obj.name = "Boundary"
        obj.index = 1
        obj.color = (255, 255, 255)

    def _fill_odata(self, bnd_types, bt, obj):
        if bt is None:
            obj.name = bnd_types.next_name()
            obj.index = bnd_types.next_index()
            obj.color = bnd_types.next_color()
        else:
            obj.name = bt.name
            obj.color = bt.color
            obj.index = bt.index

    def olist(self):
        "-> optview.OptionsList"
        return optview.OptionsList([("Boundary Type", "Index",
                optview.BoundedIntOptionEntry(self.odata(), "index",
                self._ind0, 1e6)),
            ("Boundary Type", "Name",
                optview.SimpleOptionEntry(self.odata(), "name")),
            ("Boundary Type", "Color",
                optview.ColorOptionEntry(self.odata(), "color")),
        ])

    def ret_value(self):
        '-> index, name, (char)*3 color'
        od = copy.deepcopy(self.odata())
        return od.index, od.name, od.color

    def _active_entries(self, entry):
        """ (optview.OptionEntry entry) -> bool
            return False for non-active entries
        """
        if entry.member_name in ["name", "index"]:
            return self._ind0 != 0
        else:
            return True


class GridViewOpt(QtGui.QDialog, qtui.ui_GridViewOpt.Ui_Dialog):
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
