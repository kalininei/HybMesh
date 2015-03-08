#!/usr/bin/env python
'Commands history widget'

import functools
from PyQt4 import QtCore, QtGui
import qtui.ComGridRes_rc
qtui.ComGridRes_rc


#===================== Command table MVC
class CommandTableConfigure(object):
    "Configuration for command history window"
    def __init__(self):
        #font
        self.font = QtGui.QFont()
        self.font.setPointSize(12)
        #colors
        self.bgcolor1 = QtGui.QColor(255, 255, 180)
        self.bgcolor2 = QtGui.QColor(255, 255, 200)
        self.inactive_bgcolor = QtGui.QColor(252, 252, 252)
        self.active_font_color = QtGui.QColor(0, 0, 0)
        self.inactive_font_color = QtGui.QColor(150, 150, 150)
        self.selectcolor = QtGui.QColor(38, 169, 149)
        #maximum command column width
        self.max_name_width = 200
        #maximum row height
        self.max_height = 10000
        #user comments editor size
        self.editor_size = (400, 400)
        #cellsize = textsize + text_margin
        self.text_margin = 10

comtab_opt = CommandTableConfigure()


class CommandTableDelegate(QtGui.QAbstractItemDelegate):
    """ delegate for command table:
        command name | user comments table.
        Comments are editted in QTextEdit widget
    """
    def __init__(self, header, parent=None):
        """ (QHeaderView header, QWidget parent)
            header is a horizontal header. It is used to calculate
            column widths with setColumnWidthToContent feature"""
        super(CommandTableDelegate, self).__init__(parent)
        self.header = header

    def _is_enabled(self, option):
        "extract enabled state from QStyleOptionViewItem"
        return bool(QtGui.QStyle.State_Enabled & option.state)

    class Editor(QtGui.QDialog):
        "user comments editor"
        def __init__(self, close_signal, commit_signal, parent=None):
            super(CommandTableDelegate.Editor, self).__init__(parent)
            self.close_signal = close_signal
            self.commit_signal = commit_signal
            self.tw = QtGui.QTextEdit(self)
            self.setFocusProxy(self.tw)
            bbox = QtGui.QDialogButtonBox(self)
            bbox.setStandardButtons(QtGui.QDialogButtonBox.Cancel |
                    QtGui.QDialogButtonBox.Ok)
            bbox.accepted.connect(self.accept)
            bbox.rejected.connect(self.reject)
            layout = QtGui.QVBoxLayout(self)
            layout.addWidget(self.tw)
            layout.addWidget(bbox)

        def set_text(self, txt):
            "set text to editor"
            self.tw.setText(txt)

        def get_text(self):
            "returns string from editor"
            return self.tw.toPlainText()

        def accept(self):
            "overriden"
            self.commit_signal.emit(self)
            self.close_signal.emit(self, QtGui.QAbstractItemDelegate.NoHint)
            super(CommandTableDelegate.Editor, self).accept()

        def reject(self):
            "overriden"
            self.close_signal.emit(self, QtGui.QAbstractItemDelegate.NoHint)
            super(CommandTableDelegate.Editor, self).reject()

    def createEditor(self, parent, option, index):
        "overriden"
        return self.Editor(self.closeEditor, self.commitData, parent)

    def paint(self, painter, option, index):
        "overriden"
        wdg = QtGui.QTextEdit(self.parent())
        wdg.setText(index.data(QtCore.Qt.DisplayRole).toString())
        wdg.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        wdg.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        wdg.setAutoFillBackground(True)
        wdg.setFrameShape(QtGui.QFrame.NoFrame)
        if self._is_enabled(option):
            tc = comtab_opt.active_font_color
            bc = comtab_opt.bgcolor1 if index.row() % 2 == 0 \
                    else comtab_opt.bgcolor2
        else:
            tc = comtab_opt.inactive_font_color
            bc = comtab_opt.inactive_bgcolor
        ssh1 = "color: rgba(%i, %i, %i, %i)" % tc.getRgb()
        ssh2 = "background-color: rgba(%i, %i, %i, %i)" % bc.getRgb()
        wdg.setStyleSheet("; ".join([ssh1, ssh2]))

        painter.save()
        painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
        painter.translate(option.rect.topLeft())
        wdg.resize(option.rect.size())
        wdg.render(painter)
        painter.restore()

    def setEditorData(self, editor, index):
        "overriden"
        editor.set_text(index.data(QtCore.Qt.EditRole).toString())

    def setModelData(self, editor, model, index):
        "overriden"
        t = QtCore.QVariant(editor.get_text())
        model.setData(index, t, QtCore.Qt.EditRole)

    def _widths(self, option, index):
        "->(names width, comments width). Calculate table column widths"
        metrics = option.fontMetrics
        maxr = QtCore.QRect(0, 0, comtab_opt.max_name_width,
                comtab_opt.max_height)
        lw = 0
        for i in range(index.model().rowCount()):
            ind = index.sibling(i, 0)
            t = ind.data().toString()
            brect = metrics.boundingRect(maxr, QtCore.Qt.TextWordWrap, t)
            if brect.width() > lw:
                lw = brect.width()
        lw += comtab_opt.text_margin
        rw = self.header.size().width() - lw
        return lw, rw

    def sizeHint(self, option, index):
        "overriden"
        metrics = option.fontMetrics
        lw, rw = self._widths(option, index)
        t = index.data().toString()
        if index.column() == 0:
            maxr = QtCore.QRect(0, 0, lw, comtab_opt.max_height)
            brect = metrics.boundingRect(maxr, QtCore.Qt.TextWordWrap, t)
            lh = min(brect.height(), comtab_opt.max_height)
            return QtCore.QSize(lw, lh + comtab_opt.text_margin)
        else:
            maxr = QtCore.QRect(0, 0, rw, comtab_opt.max_height)
            brect = metrics.boundingRect(maxr, QtCore.Qt.TextWordWrap, t)
            return QtCore.QSize(rw, brect.height() + comtab_opt.text_margin)

    def updateEditorGeometry(self, editor, option, index):
        "overriden"
        r = QtGui.QApplication.desktop().screenGeometry()
        x = (r.left() + r.right() - editor.geometry().width()) / 2
        y = (r.bottom() + r.top() - editor.geometry().height()) / 2
        cnt = QtCore.QPoint(x, y)
        editor.move(cnt)
        editor.resize(*comtab_opt.editor_size)


class CommandTableModel(QtCore.QAbstractTableModel):
    "command table model uses data from command.COmmandFlow"
    def __init__(self, comflow, parent=None):
        "(command.CommandFlow, QWidget parent)"
        super(CommandTableModel, self).__init__(parent)
        self.doicon = QtGui.QPixmap(":/icons/do.png")
        self.flow = None
        self.reset_comflow(comflow)

    def reset_comflow(self, comflow):
        "sets command flow to current model"
        if self.flow is not None:
            self.flow.remove_subscriber(self)
        self.flow = comflow
        self.flow.add_subscriber(self)
        self.reset()

    def command_flow_action(self, tp):
        "Executed by self.flow messager"
        self.reset()

    def rowCount(self, parent=None):
        "overriden"
        return self.flow.com_count()

    def columnCount(self, parent=None):
        "overriden"
        return 2

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        "overriden"
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Vertical:
                if section == self.flow._curpos:
                    s = ""
                else:
                    s = " %s " % str(section)
            else:
                s = " Command " if section == 0 else " Comments "
            return QtCore.QVariant(s)
        elif role == QtCore.Qt.DecorationRole:
            if orientation == QtCore.Qt.Vertical:
                if section == self.flow._curpos:
                    return self.doicon

    def data(self, index, role):
        "overriden"
        if role in [QtCore.Qt.DisplayRole, QtCore.Qt.EditRole]:
            if index.column() == 0:
                return self.flow.com(index.row()).doc()
            elif index.column() == 1:
                return self.flow.com(index.row()).get_comment()
        return None

    def setData(self, index, value, role):
        "overriden"
        if role == QtCore.Qt.EditRole:
            v = str(value.toString())
            self.flow.com(index.row()).set_comment(v)
            self.dataChanged.emit(index, index)
            return True
        return False

    def flags(self, index):
        "overriden"
        if self.flow._startpos <= index.row() <= self.flow._curpos:
            ret = QtCore.Qt.ItemIsEnabled
            if index.column() == 1:
                ret |= QtCore.Qt.ItemIsEditable
            return ret
        else:
            return QtCore.Qt.NoItemFlags


class CommandTable(QtGui.QTableView):
    "Table widget for command flow representation"
    def __init__(self, comflow, parent=None):
        "CommandTable(command.CommandFlow comflow)"
        super(CommandTable, self).__init__(parent)
        self.setVerticalScrollMode(self.ScrollPerPixel)

        hh = self.horizontalHeader()
        hh.setResizeMode(QtGui.QHeaderView.ResizeToContents)
        hh.setStretchLastSection(True)

        vh = self.verticalHeader()
        vh.setResizeMode(QtGui.QHeaderView.ResizeToContents)
        #vh.setDefaultAlignment(QtCore.Qt.AlignTop)

        delegate = CommandTableDelegate(hh, self)
        self.setItemDelegate(delegate)

        model = CommandTableModel(comflow, self)
        self.setModel(model)

    def reset_comflow(self, comflow):
        "set new CommandFlow to table"
        self.model().reset_comflow(comflow)


class FlowListModel(QtCore.QAbstractTableModel):
    """ Model for COmmands Flows representation.
        Flows are shown in a list with additional
        icon buttons: remove flow, checkpoint flow
        and context menu for each flow"""

    def __init__(self, flow_collection, parent=None):
        "(command.FlowCollection flow_collection, QWidget parent)"
        super(FlowListModel, self).__init__(parent)
        self.flows = flow_collection
        self.rem_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/delete.png"))
        self.chp_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/split.png"))

    def rowCount(self, parent=None):
        "overriden"
        return len(self.flows.get_flow_names())

    def columnCount(self, parent=None):
        "overriden"
        return 3

    def data(self, index, role):
        "overriden"
        if role in [QtCore.Qt.DisplayRole, QtCore.Qt.EditRole]:
            if index.column() == 0:
                s = self.flows.get_flow_names()[index.row()]
                #force column[0] to have at least 15 chars
                if len(s) < 15:
                    s += " " * (15 - len(s))
                return QtCore.QVariant(s)
        elif role == QtCore.Qt.DecorationRole:
            #icons for remove and checkpoint
            if index.column() == 1:
                return self.chp_icon
            if index.column() == 2:
                return self.rem_icon
        else:
            #active row is highlighted, bold font
            if index.column() == 0:
                s1 = self.flows.get_flow_names()[index.row()]
                s2 = self.flows.get_actual_flow_name()
                if s1 == s2:
                    if role == QtCore.Qt.BackgroundColorRole:
                        return comtab_opt.selectcolor
                    elif role == QtCore.Qt.FontRole:
                        font = QtGui.QFont()
                        font.setBold(True)
                        return font

    def flags(self, index):
        "overriden"
        return QtCore.Qt.ItemIsEnabled


class FlowList(QtGui.QTableView):
    "Widget for flow representation"
    def __init__(self, flow_collection, parent=None):
        "(command.FlowCollection flow_collection, QWidget parent)"
        super(FlowList, self).__init__(parent)
        self.flow_collection = flow_collection
        self.model = FlowListModel(flow_collection, self)
        self.setModel(self.model)

        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setResizeMode(
                QtGui.QHeaderView.ResizeToContents)
        self.verticalHeader().setVisible(False)

        self.doubleClicked.connect(self._double_click)
        self.clicked.connect(self._click)

        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self._context_menu)

        self.setShowGrid(False)

    def _context_menu(self, pnt):
        index = self.indexAt(pnt)
        if index.isValid() and index.column() == 0:
            def trig(m, name, fun, i):
                "create menu item"
                act = QtGui.QAction(name, self)
                act.triggered.connect(functools.partial(fun, i))
                m.addAction(act)
            menu = QtGui.QMenu(self)
            trig(menu, "Set active", self._set_active, index.row())
            trig(menu, "Rename", self._rename, index.row())
            trig(menu, "Checkpoint", self._checkpoint, index.row())
            trig(menu, "To start", self._to_start, index.row())
            trig(menu, "Remove", self._remove, index.row())
            menu.popup(self.viewport().mapToGlobal(pnt))

    def _count(self):
        return len(self.flow_collection.get_flow_names())

    def _double_click(self, index):
        if index.column() == 0:
            self._set_active(index.row())

    def _click(self, index):
        if index.column() == 0:
            return
        if index.column() == 2:
            self._remove(index.row())
        elif index.column() == 1:
            self._checkpoint(index.row())

    def _set_active(self, i):
        fn = self.flow_collection.get_flow_names()[i]
        self.flow_collection.set_actual_flow(fn)

    def _remove(self, i):
        fn = self.flow_collection.get_flow_names()[i]
        a = QtGui.QMessageBox.question(None, "Confirmation",
                "Delete flow " + fn + "?",
                QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if a == QtGui.QMessageBox.Yes:
            try:
                self.flow_collection.remove_flow(fn)
                self.reset()
                self.resizeColumnsToContents()
            except:
                QtGui.QMessageBox.warning(None, "Warning",
                        "Unable to remove %s" % fn)

    def _to_start(self, i):
        self.flow_collection.get_flow(ind=i).to_zero_state()

    def _rename(self, i):
        fn = self.flow_collection.get_flow_names()[i]
        nm, ok = QtGui.QInputDialog.getText(self, "Rename",
                "Enter new flow name", text=fn)
        fn, nm = str(fn), str(nm)
        if ok and fn != nm:
            self.flow_collection.rename_flow(fn, str(nm))
            self.reset()
            self.resizeColumnsToContents()

    def _checkpoint(self, i):
        fn = self.flow_collection.get_flow_names()[i]
        nm, ok = QtGui.QInputDialog.getText(self, "New Checkpoint",
                "Enter checkpoint name", text="Copy_" + fn)
        if ok:
            self.model.beginInsertRows(QtCore.QModelIndex(),
                    self._count(), self._count())
            self.flow_collection.checkpoint_current(str(nm))
            self.model.endInsertRows()
            self.resizeColumnsToContents()


class HistoryDockFrame(QtGui.QFrame):
    "The Widget for commands history representation"
    def __init__(self, flow_collection, parent):
        "(command.FlowCollection flow_collection, QWidget parent)"
        super(HistoryDockFrame, self).__init__(parent)
        self.flow_collection = flow_collection
        flow_collection.add_actflow_subscriber(self)
        layout = QtGui.QHBoxLayout(self)
        splitter = QtGui.QSplitter(QtCore.Qt.Horizontal, self)
        splitter.setChildrenCollapsible(False)
        layout.addWidget(splitter)
        self.flist = FlowList(flow_collection)

        self.comtab = CommandTable(
            flow_collection.get_actual_flow(), self)
        splitter.addWidget(self.flist)
        splitter.addWidget(self.comtab)

        splitter.setSizes([150])
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

    def actual_flow_changed(self, newflow):
        "messager from FlowCollection subscription"
        #TODO: flist and comtab models should be
        #      set as subscribers?
        #list change
        self.flist.model.reset()
        #table change
        self.comtab.reset_comflow(self.flow_collection.get_actual_flow())

if __name__ == "__main__":
    import sys

    class DummyCom:
        def __init__(self, doc, com=""):
            self._com = com
            self._doc = doc

        def __str__(self):
            return "Command {arguments}"

        def doc(self):
            return self._doc

        def get_comment(self):
            return self._com

        def set_comment(self, com):
            self._com = com

    class DFlow(object):
        def add_subscriber(self, obj):
            pass

        def remove_subscriber(self, obj):
            pass

        def com_count(self):
            return len(self._commands)

        def com(self, i):
            return self._commands[i]

    class DummyFlow(DFlow):
        def __init__(self):
            self._commands = [DummyCom("Start"),
                    DummyCom("Com1", "SomeCom\n\nSomeComm3"),
                    DummyCom("Com2", "Commentary\n" * 4),
                    DummyCom("Com3", "I ad fa ds asdf asd fas ")]
            self._startpos = 1
            self._curpos = 2

    class DummyFlow2(DFlow):
        def __init__(self):
            self._commands = [DummyCom("Start"),
                    DummyCom("Command", "No Comments"),
                    DummyCom("KOMMAND", "NASTYA\n" * 4),
                    DummyCom("CAMANTA", "Ich liebe dich")]
            self._startpos = 2
            self._curpos = 3

    class DummyFlow3(DFlow):
        def __init__(self):
            self._commands = [DummyCom("Start")]
            self._startpos = 0
            self._curpos = 0

    class DummyFlowCollection:
        def __init__(self):
            self.flows = ["Flow1", "Flow2", "Flow3"]
            self.act = 0
            self.actflow_subscribers = []
            self.c = [DummyFlow2(), DummyFlow(), DummyFlow3()]

        def add_actflow_subscriber(self, obj):
            self.actflow_subscribers.append(obj)

        def get_actual_flow(self):
            return self.c[self.act]

        def get_flow_names(self):
            return self.flows[:]

        def rename_flow(self, oldname, newname):
            i = self.flows.index(oldname)
            self.flows[i] = newname

        def remove_flow(self, name):
            self.flows.remove(name)

        def checkpoint_current(self, name):
            self.flows.append(name)

        def get_actual_flow_name(self):
            return self.flows[self.act]

        def set_actual_flow(self, flow_name):
            self.act = self.flows.index(flow_name)
            for o in self.actflow_subscribers:
                o.actual_flow_changed(flow_name)

    class Window(QtGui.QWidget):
        def __init__(self, parent=None):
            super(Window, self).__init__(parent)
            tab = HistoryDockFrame(DummyFlowCollection(), self)
            layout = QtGui.QVBoxLayout(self)
            layout.addWidget(tab)

    app = QtGui.QApplication(sys.argv)
    win = Window()
    win.show()
    sys.exit(app.exec_())
