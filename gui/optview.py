#!/usr/bin/env python
""" Tree like options representation widget.

Defined public classes:
1) Basic classes
    OptionsView - widget derived from QTreeView

    OptionsList - driver which sets a match between
        value member in a user class and option caption/name and display/edit
        representation in the widget

    OptionWidgetConfigure - color, font, indent widget configuation.
        Default configuration is stored in OptionWidgetConfigure.default.

2) Option representation classes (are used for OptionsList building):
    SimpleOptionEntry - QLabel (for display) + QLineEdit (for edit).
        Can be used for arbitrary str, int, float data.

    BoundedIntOptionEntry - QLabel + QSpinEdit.
        Used for integers with defined bounds

    BoolOptionEntry - QCheckBox + QCheckBox. Boolean entries.

    SingleChoiceOptionEntry - QLabel + QComboBox.
        Pick a single string from a string set

    MultipleChoiseOptionEntry - Multiple QLabels + External widget
        Pick sublist of strings from a string set

    XYOptionEntry - QLabel + pair of QLineEdits.
        Used for classes with x, y properties
"""

import sys
import copy
import collections
from PyQt4 import QtCore, QtGui
import qtui.ui_ComGridMain
qtui.ui_ComGridMain


# ================ Configurations
class OptionWidgetConfigure(object):
    "defines color, font etc. of option table widget "
    def __init__(self):
        #colors
        self.bgcolor1 = QtGui.QColor(255, 255, 0, 70)
        self.bgcolor2 = QtGui.QColor(255, 255, 0, 40)
        self.captioncolor = QtGui.QColor(150, 150, 150, 100)
        self.selectcolor = QtGui.QColor(38, 169, 149)
        self.active_font_color = QtGui.QColor(0, 0, 0)
        self.inactive_font_color = QtGui.QColor(200, 200, 200)
        #font
        self.font = QtGui.QFont()
        self.font.setPointSize(10)
        #cell standard height
        self.row_height = 20
        #indentation
        self.font_indent = 5
        self.branch_indent = 20
        self.vert_gap = 2

    def _string_width(self, s):
        m = QtGui.QFontMetrics(self.font)
        return m.width(s)

    def _label(self, text="", parent=None):
        wdg = QtGui.QLabel(text, parent)
        wdg.setIndent(self.font_indent)
        wdg.setFont(self.font)
        return wdg

OptionWidgetConfigure.default = OptionWidgetConfigure()


# ================ Option Entries
class OptionEntry(object):
    "Options entry base class "
    def __init__(self, data, member_name):
        """ (OptionsSet data, String dict_name) """
        self.data = data
        self.member_name = member_name
        self._check(self.get())
        self.conf = OptionWidgetConfigure.default

    def set_configuration(self, conf):
        "set OptionWidgetConfigure for the item"
        self.conf = conf

    def _check(self, val):
        "throws ValueError if val is not correct"
        if not self._check_proc(val):
            raise ValueError

    def _check_proc(self, val):
        "function for overriding: -> True if value is ok"
        return True

    #widgets
    def display_widget(self):
        " generates display widget with data "
        raise NotImplementedError

    def edit_widget(self, parent=None):
        " genenerates edit widget with initial data"
        raise NotImplementedError

    def display_lines(self):
        "number of lines in display widget is used to define cell size"
        return 1

    def edit_lines(self):
        "number of lines in edit widget is used to define cell size"
        return 1

    #values manipulations
    def get(self):
        return self.data.__dict__[self.member_name]

    def set(self, val):
        'check value and pass it to data'
        self._check(val)
        self.data.__dict__[self.member_name] = val

    def set_from_widget(self, editor):
        "set value from edit_widget"
        raise NotImplementedError

    #value to string
    def __str__(self):
        return str(self.get())


#option entries types
class SimpleOptionEntry(OptionEntry):
    "simple str, int, float option entries "
    def __init__(self, data, member_name):
        super(SimpleOptionEntry, self).__init__(data, member_name)
        self.tp = type(self.get())

    def display_widget(self):
        return self.conf._label(str(self))

    def edit_widget(self, parent):
        wdg = QtGui.QLineEdit(parent)
        wdg.setText(str(self))
        return wdg

    def set_from_widget(self, widget):
        self.set(self.tp(widget.text()))


class BoundedIntOptionEntry(SimpleOptionEntry):
    " integer value within [minv, maxv] "
    def __init__(self, data, member_name, minv, maxv):
        self.minv = minv
        self.maxv = maxv
        super(BoundedIntOptionEntry, self).__init__(data, member_name)

    def _check_proc(self, v):
        return v >= self.minv and v <= self.maxv

    def edit_widget(self, parent):
        wdg = QtGui.QSpinBox(parent)
        wdg.setMinimum(self.minv)
        wdg.setMaximum(self.maxv)
        wdg.setValue(self.get())
        return wdg

    def set_from_widget(self, widget):
        self.set(widget.value())


class SingleChoiceOptionEntry(SimpleOptionEntry):
    "string value from combobox"
    def __init__(self, data, member_name, values):
        self.values = copy.deepcopy(values)
        super(SingleChoiceOptionEntry, self).__init__(data, member_name)

    def _check_proc(self, v):
        return v in self.values

    def edit_widget(self, parent):
        wdg = QtGui.QComboBox(parent)
        wdg.addItems(self.values)
        wdg.setCurrentIndex(self.values.index(self.get()))
        return wdg

    def set_from_widget(self, widget):
        self.set(widget.currentText())


class BoolOptionEntry(SimpleOptionEntry):
    " boolean flag option"
    def __init__(self, data, member_name):
        super(BoolOptionEntry, self).__init__(data, member_name)

    def _chbox_widget(self, parent=None):
        wdg = QtGui.QCheckBox(parent)
        wdg.setChecked(self.get())
        return wdg

    def display_widget(self):
        return self._chbox_widget()

    def edit_widget(self, parent):
        return self._chbox_widget(parent)

    def set_from_widget(self, widget):
        self.set(widget.isChecked())


class XYOptionEntry(SimpleOptionEntry):
    "x, y point option entry"
    def __init__(self, data, member_name):
        super(XYOptionEntry, self).__init__(data, member_name)
        self.__last_paint = True

    def _check_proc(self, val):
        try:
            return isinstance(val.x, float) and isinstance(val.y, float)
        except:
            return False

    def __str__(self):
        v = self.get()
        return "[" + str(v.x) + ",  " + str(v.y) + "]"

    def display_lines(self):
        #1 in display mode and 3 in edit mode
        return 1 if self.__last_paint else 3

    class EditWidget(QtGui.QWidget):
        " pair of QLineEdits with X, Y labels"
        def __init__(self, x, y, parent):
            super(XYOptionEntry.EditWidget, self).__init__(parent)
            self.xed = QtGui.QLineEdit(str(x))
            self.yed = QtGui.QLineEdit(str(y))
            self.setFocusProxy(self.xed)
            layout = QtGui.QGridLayout()
            layout.addWidget(QtGui.QLabel("X"), 0, 0)
            layout.addWidget(self.xed, 0, 1)
            layout.addWidget(QtGui.QLabel("Y"), 1, 0)
            layout.addWidget(self.yed, 1, 1)
            layout.setVerticalSpacing(0)
            self.setLayout(layout)
            self.setAutoFillBackground(True)

    def display_widget(self):
        self.__last_paint = True
        return super(XYOptionEntry, self).display_widget()

    def edit_widget(self, parent):
        self.__last_paint = False
        v = self.get()
        return XYOptionEntry.EditWidget(v.x, v.y, parent)

    def set_from_widget(self, widget):
        try:
            x = float(widget.xed.text())
            y = float(widget.yed.text())
            self.get().x = x
            self.get().y = y
        except:
            QtCore.qDebug("Invalid x, y value")


class MultipleChoiceOptionEntry(OptionEntry):
    'Unique sublist from a list of string. Ordering matters.'

    def __init__(self, data, member_name, values):
        self.values = copy.deepcopy(values)
        super(MultipleChoiceOptionEntry, self).__init__(data, member_name)

    def _check_proc(self, val):
        return all(v in self.values for v in val)

    def display_widget(self):
        " generates display widget with data "
        wdg = QtGui.QFrame()
        lab = QtGui.QVBoxLayout()
        lab.setSpacing(self.conf.vert_gap)
        lab.setMargin(0)
        lab.addStretch(1)
        for v in self.get():
            lab.addWidget(self.conf._label(v))
        lab.addStretch(1)
        wdg.setLayout(lab)
        return wdg

    def display_lines(self):
        return len(self.get())

    class EditWidget(QtGui.QWidget):
        'widget for creating unique sublist'
        def __init__(self, data, parent):
            #---- building window
            super(MultipleChoiceOptionEntry.EditWidget, self).__init__(parent)
            self.lw1 = QtGui.QListWidget()
            self.lw2 = QtGui.QListWidget()
            btleft = QtGui.QPushButton()
            btleft.setIcon(QtGui.QIcon(
                QtGui.QPixmap(":/icons/left_arrow.png")))
            btleft.clicked.connect(self.left_click)
            btleft.setFocusPolicy(QtCore.Qt.NoFocus)
            btright = QtGui.QPushButton()
            btright.clicked.connect(self.right_click)
            btright.setIcon(QtGui.QIcon(
                QtGui.QPixmap(":/icons/right_arrow.png")))
            btright.setFocusPolicy(QtCore.Qt.NoFocus)
            bbox = QtGui.QDialogButtonBox()
            bbox.setStandardButtons(QtGui.QDialogButtonBox.Cancel |
                    QtGui.QDialogButtonBox.Ok)
            bbox.accepted.connect(self.ok_action)
            bbox.rejected.connect(self.cancel_action)
            #left/right buttons
            button_frame = QtGui.QFrame()
            button_frame.setLayout(QtGui.QVBoxLayout())
            button_frame.layout().addStretch(1.0)
            button_frame.layout().addWidget(btleft)
            button_frame.layout().addWidget(btright)
            button_frame.layout().addStretch(1.0)
            #widget
            widget_frame = QtGui.QFrame()
            widget_frame.setLayout(QtGui.QHBoxLayout())
            widget_frame.layout().addWidget(self.lw1)
            widget_frame.layout().addWidget(button_frame)
            widget_frame.layout().addWidget(self.lw2)
            #window
            self.setFocusProxy(self.lw1)
            self.setLayout(QtGui.QVBoxLayout())
            self.layout().addWidget(widget_frame)
            self.layout().addWidget(bbox)
            self.resize(400, 300)
            self.setWindowModality(QtCore.Qt.WindowModal)
            #---- data
            self.data = data
            self.left_column = data.get()
            self.right_column = filter(lambda x: x not in self.left_column,
                    data.values)
            self._fill()

        def _fill_column(self, lw, col, s):
            lw.clear()
            map(lw.addItem, col)
            if len(col) > 0:
                s = min(s, len(col) - 1)
                lw.setCurrentRow(s)

        def _fill(self, s1=0, s2=0):
            """fill widget lists from data arrays
                s1, s2 - selected row indicies
            """
            self._fill_column(self.lw1, self.left_column, s1)
            self._fill_column(self.lw2, self.right_column, s2)

        def _add_rem(self, remcol, addcol, itms):
            s1, s2 = self.lw1.currentRow(), self.lw2.currentRow()
            for it in [str(v.text()) for v in itms]:
                remcol.remove(it)
                addcol.append(it)
            self._fill(s1, s2)

        def left_click(self):
            self._add_rem(self.right_column, self.left_column,
                    self.lw2.selectedItems())

        def right_click(self):
            self._add_rem(self.left_column, self.right_column,
                    self.lw1.selectedItems())

        def ok_action(self):
            self.data.set(self.left_column)
            self.close()

        def cancel_action(self):
            self.close()

        def keyPressEvent(self, event):
            k = event.key()
            if self.focusWidget() is self.lw1 and k == QtCore.Qt.Key_Right:
                return self.right_click()
            elif self.focusWidget() is self.lw2 and k == QtCore.Qt.Key_Left:
                return self.left_click()
            elif k in [QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter]:
                return self.ok_action()
            super(self.__class__, self).keyPressEvent(event)

    def edit_widget(self, parent=None):
        #create widget in a separate window without parents
        return MultipleChoiceOptionEntry.EditWidget(self, None)

    def get(self):
        v = super(MultipleChoiceOptionEntry, self).get()
        return copy.deepcopy(v)

    def set(self, val):
        v = copy.deepcopy(val)
        super(MultipleChoiceOptionEntry, self).set(v)

    def set_from_widget(self, editor):
        #value was set in editor widget ok action
        pass

    def __str__(self):
        return "; ".join(self.get())


class OptionsList(object):
    'list of objects derived from OptionEntry class'

    def __init__(self, opt_array):
        """ should be initialized from the list of tuples:
            ([(Caption1, OptionName1, OptionEntry1),
            (Caption2, OptionName2, OptionEntry2),...]
        """
        self.caps = [x[0] for x in opt_array]
        self.names = [x[1] for x in opt_array]
        self.opts = [x[2] for x in opt_array]

    def set_configuration(self, conf):
        'set OptionWidgetConfigure for all entries'
        for v in self.opts:
            v.set_configuration(conf)

    def captions(self):
        " -> list of captions"
        return list(collections.OrderedDict.fromkeys(self.caps))

    def cap_options(self, caption):
        """-> [(OptionName1, OptionEntry1), (OptionName2, OptionEntry2), ...]
        returns set of options for certain caption
        """
        return [(y, z) for x, y, z in zip(self.caps, self.names, self.opts)
                    if x == caption]

    def __str__(self):
        lines = ["------- OPTIONS"]
        for i in range(len(self.opts)):
            lines.append(self.names[i] + ": " + str(self.opts[i]))
        return '\n'.join(lines)


# ============================ Delegate
class OptionsValueDelegate(QtGui.QItemDelegate):
    'representation of options in a widget'
    def __init__(self, conf):
        super(OptionsValueDelegate, self).__init__()
        self.conf = conf

    #--- index categories
    def _index_position(self, index):
        """-> int:
            0 - no data, 1 - group caption,
            2 - option name, 3 - option data
        """
        if index.column() == 0:
            return 2 if index.internalPointer().childCount() == 0 else 1
        else:
            return 3 if index.internalPointer().data(1) is not None else 0

    def _bg_color(self, index):
        "bg color for the tree row"
        if self._index_position(index) in [0, 1]:
            return self.conf.captioncolor
        else:
            return self.conf.bgcolor1 if index.row() % 2 == 0 \
                    else self.conf.bgcolor2

    def _option_data(self, index):
        "option value from the index"
        return index.internalPointer().data(1)

    #--- options categories
    def _is_checked(self, option):
        "extract check state from QStyleOptionViewItem"
        return bool(QtGui.QStyle.State_Selected & option.state)

    def _is_enabled(self, option):
        "extract enabled state from QStyleOptionViewItem"
        return bool(QtGui.QStyle.State_Enabled & option.state)

    def _caption_widget(self, index):
        "build a label widget for the index row"
        pos = self._index_position(index)
        wdg = self.conf._label(index.data().toString())
        f = QtGui.QFont(self.conf.font)
        f.setBold(pos == 1)
        wdg.setFont(f)
        return wdg

    def _widget_color(self, option):
        "-> QColor. transparent or selection color"
        if self._is_checked(option):
            return self.conf.selectcolor
        else:
            return QtGui.QColor(0, 0, 0, 0)

    def _widget_font_color(self, option):
        "-> QColor. font color depending on is_enabled property"
        if self._is_enabled(option):
            return self.conf.active_font_color
        else:
            return self.conf.inactive_font_color

    def paint(self, painter, option, index):
        "overriden"
        pos = self._index_position(index)
        #gridlines and color
        pen = QtGui.QPen()
        pen.setColor(QtCore.Qt.lightGray)
        painter.setPen(pen)
        rect = QtCore.QRect(option.rect)
        if pos in [1, 2]:
            rect.setLeft(0)
        painter.fillRect(rect, QtGui.QBrush(self._bg_color(index)))
        if pos in [2, 3]:
            painter.drawRect(rect)
        else:
            painter.drawLine(rect.topLeft(), rect.topRight())

        #display widgets
        if pos == 0:
            #blank square
            super(OptionsValueDelegate, self).paint(painter, option, index)
        else:
            #construct widget
            if pos in [1, 2]:
                #group/option caption
                wdg = self._caption_widget(index)
            else:
                #option value
                wdg = self._option_data(index).display_widget()

            #set transparency or checked color for widget
            #set font color enabled or disabled
            ssh1 = "background-color: rgba(%i, %i, %i, %i)" \
                    % self._widget_color(option).getRgb()
            ssh2 = "color: rgba(%i, %i, %i, %i)"  \
                    % self._widget_font_color(option).getRgb()
            wdg.setStyleSheet(";".join([ssh1, ssh2]))
            #paint widget
            painter.save()
            painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
            painter.translate(option.rect.topLeft())
            wdg.resize(option.rect.size())
            wdg.render(painter)
            painter.restore()

    def createEditor(self, parent, option, index):
        "overriden"
        ed = self._option_data(index).edit_widget(parent)
        #emit size hint in case editor has more lines then display widget
        self.sizeHintChanged.emit(index)
        return ed

    def sizeHint(self, option, index):
        "overriden"
        pos = self._index_position(index)
        if pos == 3:
            addl = max(1, self._option_data(index).display_lines()) - 1
            delta = (self.conf.row_height - 2 * self.conf.vert_gap) * addl
            return QtCore.QSize(-1, self.conf.row_height + delta)
        elif pos in [1, 2]:
            w = self.conf._string_width(index.data().toString())
            w += 2 * self.conf.branch_indent + self.conf.font_indent
            return QtCore.QSize(w, self.conf.row_height)
        else:
            return QtCore.QSize(-1, -1)

    def setModelData(self, editor, model, index):
        "overriden"
        self._option_data(index).set_from_widget(editor)
        #emit size hint in case of change of display lines
        self.sizeHintChanged.emit(index)

    def updateEditorGeometry(self, editor, option, index):
        "overriden"
        #place external editor in the center of the screen
        if editor.parent() is None:
            r = QtGui.QApplication.desktop().screenGeometry()
            x = (r.left() + r.right() - editor.geometry().width()) / 2
            y = (r.bottom() + r.top() - editor.geometry().height()) / 2
            cnt = QtCore.QPoint(x, y)
            editor.move(cnt)
        else:
            super(OptionsValueDelegate, self).updateEditorGeometry(editor,
                    option, index)


# ============================= Model
class TreeItem(object):
    """item in a option tree. Base abstract class.
    Item is represented as a row in a QTreeView"""

    def __init__(self, data, parent=None):
        self.parentItem = parent
        self.itemData = data
        self.childItems = []
        if parent is not None:
            parent.appendChild(self)

    def appendChild(self, item):
        self.childItems.append(item)

    def child(self, row):
        return self.childItems[row]

    def childCount(self):
        return len(self.childItems)

    def columnCount(self):
        return 2

    def data(self, column):
        return None

    def parent(self):
        return self.parentItem

    def row(self):
        return self.parentItem.childItems.index(self) \
                if self.parentItem else 0


class CaptionTreeItem(TreeItem):
    "Caption row in a treeview"
    def __init__(self, data, parent=None):
        "(caption string, root item)"
        super(CaptionTreeItem, self).__init__(data, parent)

    def data(self, column):
        return self.itemData if column == 0 else None


class ValueTreeItem(TreeItem):
    "name: option row in a tree view"
    def __init__(self, data, parent):
        "((name, option entry), CaptionTreeItem)"
        super(ValueTreeItem, self).__init__(data, parent)

    def data(self, column):
        return self.itemData[column]


class OptionsModel(QtCore.QAbstractItemModel):
    "Tree like model for options representation"
    def __init__(self, opts):
        super(OptionsModel, self).__init__()
        self._root_item = TreeItem("Root")
        self._setup_model_data(opts, self._root_item)
        self.is_active = lambda x: True

    def rowCount(self, parent):
        "overriden"
        if parent.column() > 0:
            return 0

        if not parent.isValid():
            return self._root_item.childCount()
        else:
            return parent.internalPointer().childCount()

    def columnCount(self, parent):
        "overriden"
        if parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return self._root_item.columnCount()

    def index(self, row, column, parent):
        "overriden"
        if not self.hasIndex(row, column, parent):
            return QtCore.QModelIndex()

        if not parent.isValid():
            child = self._root_item.child(row)
        else:
            child = parent.internalPointer().child(row)

        if child:
            return self.createIndex(row, column, child)
        else:
            return QtCore.QModelIndex()

    def parent(self, index):
        "overriden"
        if not index.isValid():
            return QtCore.QModelIndex()

        parent = index.internalPointer().parent()
        if parent == self._root_item:
            return QtCore.QModelIndex()

        return self.createIndex(parent.row(), 0, parent)

    def data(self, index, role):
        "overriden"
        if not index.isValid():
            return None
        item = index.internalPointer()
        if role == QtCore.Qt.DisplayRole or role == QtCore.Qt.EditRole:
            return item.data(index.column())
        else:
            return None

    def flags(self, index):
        "overriden"
        item = index.internalPointer()
        ret = QtCore.Qt.ItemIsEnabled
        #only value rows are considered
        if item.childCount() == 0:
            if self.is_active(item.data(1)):
                #active rows are selectable and editable
                ret = ret | QtCore.Qt.ItemIsSelectable
                if index.column() == 1:
                    ret = ret | QtCore.Qt.ItemIsEditable
            else:
                #no flags for disabled rows
                ret = QtCore.Qt.NoItemFlags
        return ret

    def _setup_model_data(self, data, root):
        "builds a model tree"
        caps = data.captions()
        for c in caps:
            newc = CaptionTreeItem(c, root)
            for v in data.cap_options(c):
                ValueTreeItem(v, newc)


#---------------------- OptionsView
class OptionsView(QtGui.QTreeView):
    "options representation widget"
    def __init__(self, opt_list, conf=OptionWidgetConfigure.default,
            parent=None):
        "(OptionsList opt_list, OptionWidgetConfigure conf, QWidget parent)"
        super(OptionsView, self).__init__(parent)
        #apply configuration
        opt_list.set_configuration(conf)
        #build model/view
        self.model = OptionsModel(opt_list)
        self.delegate = OptionsValueDelegate(conf)
        self.setModel(self.model)
        self.setItemDelegate(self.delegate)
        #set QTreeView options
        self.expandAll()
        self.resizeColumnToContents(0)
        self.setAllColumnsShowFocus(True)
        self.setHeaderHidden(True)
        self.setIndentation(conf.branch_indent)

    def is_active_delegate(self, func):
        """ define a function which sets active status for entries
                func = bool function(OptionEntry)
        """
        self.model.is_active = func

    def keyPressEvent(self, event):
        "overriden"
        #add enter (return key in qt notation) pressed to editor event
        if event.key() in [QtCore.Qt.Key_Enter, QtCore.Qt.Key_Return]:
            for index in self.selectedIndexes():
                if index.flags() & QtCore.Qt.ItemIsEditable:
                    trigger = QtGui.QAbstractItemView.EditKeyPressed
                    self.edit(index, trigger, event)
                    break
        else:
            #using "else:" here because otherwise it emits enter signal
            #for the whole form
            super(OptionsView, self).keyPressEvent(event)


#=================== Usage example
if __name__ == "__main__":
    class Point(object):
        def __init__(self, x, y):
            self.x, self.y = x, y

    class CircOptionsSet:
        def __init__(self):
            self.name = "CircularGrid1"
            self.px = 0.0
            self.py = 0.0
            self.rad = 1.0
            self.Nr = 10
            self.Na = 10
            self.is_trian = True
            self.units = "units"
            self.grids = ["Grid1", "Grid2", "Grid3"]
            self.point = Point(3.4, 6.7)

    def optionlist_for_circ(opt):
        #all possible values for opt.units, opt.grids
        a_units = ["units", "%"]
        a_grids = ["Grid1", "Grid2", "Grid3", "Grid4", "Grid5"]
        #build an input array for OptionsList:
        #   [(Caption, OptionName, OptionEntry), ....]
        ar = [
                ("Basic", "Grid name", SimpleOptionEntry(opt, "name")),
                ("Geometry", "X coordinate", SimpleOptionEntry(opt, "px")),
                ("Geometry", "Y coordinate", SimpleOptionEntry(opt, "py")),
                ("Geometry", "Radius", SimpleOptionEntry(opt, "rad")),
                ("Partition", "Radius Partition",
                    BoundedIntOptionEntry(opt, "Nr", 1, 1e2)),
                ("Partition", "Arch Partition",
                    BoundedIntOptionEntry(opt, "Na", 3, 1e2)),
                ("Partition", "Triangulate center cell",
                    BoolOptionEntry(opt, "is_trian")),
                ("Additional", "Units",
                    SingleChoiceOptionEntry(opt, "units", a_units)),
                ("Additional", "Choice",
                    MultipleChoiceOptionEntry(opt, "grids", a_grids)),
                ("Additional", "Point", XYOptionEntry(opt, "point"))
            ]
        return OptionsList(ar)

    class Window(QtGui.QDialog):
        def __init__(self, opt, parent=None):
            super(Window, self).__init__(parent)
            self.resize(400, 500)
            self.button = QtGui.QPushButton()
            self.button.setText("Options to stdout")
            #!Option widget
            self.opt_list = optionlist_for_circ(opt)
            self.tab = OptionsView(self.opt_list, parent=self)

            #add smart active row management
            def is_active(entry):
                #matches name entry activity with is_trian flag
                if entry.member_name == "name":
                    return entry.data.is_trian
                else:
                    return True
            self.tab.is_active_delegate(is_active)

            #window and layout
            vert_layout = QtGui.QVBoxLayout()
            vert_layout.addWidget(self.tab)
            vert_layout.addWidget(self.button)

            self.setLayout(vert_layout)

            self.button.clicked.connect(self.buttonClicked)

        def buttonClicked(self):
            print self.opt_list

    app = QtGui.QApplication(sys.argv)
    win = Window(CircOptionsSet())
    win.show()
    sys.exit(app.exec_())
