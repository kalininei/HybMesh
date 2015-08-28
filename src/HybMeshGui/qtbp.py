"basic qt procedures"

from PyQt4 import QtCore, QtGui
import qtui.ComGridRes_rc
qtui.ComGridRes_rc


_icon_set = {}


def get_icon(s):
    """ ->QtGui.QIcon. Returns icon by its string code:
        opts, del, eye-on, eye-off, split, ->, <-, hybmesh
    """
    global _icon_set
    if len(_icon_set) == 0:
        _icon_set = {
                "opts": QtGui.QIcon(QtGui.QPixmap(":/icons/opts.png")),
                "del": QtGui.QIcon(QtGui.QPixmap(":/icons/delete.png")),
                "eye-on": QtGui.QIcon(QtGui.QPixmap(":/icons/eye-on.png")),
                "eye-off": QtGui.QIcon(QtGui.QPixmap(":/icons/eye-off.png")),
                "split": QtGui.QIcon(QtGui.QPixmap(":/icons/split.png")),
                "->": QtGui.QIcon(QtGui.QPixmap(":/icons/right_arrow.png")),
                "<-": QtGui.QIcon(QtGui.QPixmap(":/icons/left_arrow.png")),
                "hybmesh": QtGui.QIcon(QtGui.QPixmap(":/icons/mainwin.png"))
        }
    return _icon_set[s]


class TreeItem(object):
    """item in a option tree. Base abstract item class
    Item is represented as a row in a RowsTreeModel"""

    def __init__(self, data, parent=None):
        """ data - tuple or list of represented data
            parent - TreeItem object

            Length of data should match model column count.
            Extend it with None if nesessary.
        """
        self.parent_item = parent
        self.item_data = data
        self.child_items = []
        if parent is not None:
            parent.append_child(self)

    def append_child(self, item):
        self.child_items.append(item)

    def clear_childs(self):
        del self.child_items[:]

    def child(self, row):
        return self.child_items[row]

    def child_count(self):
        return len(self.child_items)

    def column_count(self):
        return len(self.item_data)

    def parent(self):
        return self.parent_item

    def row(self):
        return self.parent_item.child_items.index(self) \
                if self.parent_item else 0

    def data(self, column):
        return self.item_data[column]


class RowsTreeModel(QtCore.QAbstractItemModel):
    """ Tree like model for data rows representation
        parent(), index(), columnCount(), rowCount() are implemented.

        Subclasses should implement at least data(self, index, role).
        flags(index) and setData(self, index, role) for editable models.
        Tree should be filled using qtbp.TreeItem in subclass constructor.
        Use self._root_item as the root qtbp.TreeItem
    """
    def __init__(self, col_count, parent=None):
        "col_count - number of columns, parent - QObject"
        super(RowsTreeModel, self).__init__(parent)
        self.col_count = col_count
        self._root_item = TreeItem(self._treeitem_data(None))

    def rowCount(self, parent):
        "overriden"
        if parent.column() > 0:
            return 0
        if not parent.isValid():
            return self._root_item.child_count()
        else:
            return parent.internalPointer().child_count()

    def columnCount(self, parent):
        "overriden"
        if parent.isValid():
            return parent.internalPointer().column_count()
        else:
            return self._root_item.column_count()

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

    def _treeitem_data(self, *args):
        """ (arg1, arg2, ...) -> (arg1, arg2, ..., None, None)
            Lenght of the returning tuple equls self.col_count.
            this is used for correct TreeItem construction
        """
        if len(args) > self.col_count:
            return args[:self.col_count]
        elif len(args) < self.col_count:
            return args + (None, ) * (self.col_count - len(args))
        else:
            return args
