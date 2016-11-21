"boundary types list"
import xml.etree.ElementTree as ET
import hybmeshpack.basic.proc as bp


class BndType(object):
    'name, color pair'
    def __init__(self, index, name, color):
        self.index = index
        self.name = name
        self.color = color


class BndTypesList(object):
    _st_cols = [(0, 0, 0),
                (50, 50, 50),
                (150, 150, 150),
                (10, 255, 0)
                ]

    'list of (index, name, color) entries'
    def __init__(self):
        self._data = [BndType(0, "None", (255, 255, 255))]

    def add_data(self, bt):
        'adds values from another BndTypesList object'
        for b in bt._data:
            self.set_bnd(b.index, b.name, b.color)

    def _ind_set(self):
        return set([b.index for b in self._data])

    def _name_set(self):
        return set([b.name for b in self._data])

    def _unique_name(self, name):
        return bp.unique_name(name, self._name_set())

    def next_index(self):
        '->next free index'
        st = self._ind_set()
        for i in range(1, int(1e3)):
            if i not in st:
                return i
        return 1

    def next_name(self):
        '->builds default name for next item'
        return self._unique_name("boundary1")

    def next_color(self):
        '->(char)*3. default color for next item'
        try:
            from collections import OrderedDict
            st_used = OrderedDict([(c, 0) for c in self._st_cols])
            used = [c.color for c in self._data]
            for c in used:
                if c in st_used:
                    st_used[c] += 1
            minused = min(st_used.values())
            return bp.find(lambda x: x[1] == minused, st_used.items())[0]
        except:
            return (255, 255, 255)

    def bnd_count(self):
        return len(self._data)

    def clear(self):
        self.__init__()

    def set_bnd(self, index, name, color):
        """ (int index, str name, (tuple-of-chars) color)
            if index exists -> overwrites,
            if name exists -> changes name
        """
        import operator
        if index in self._ind_set():
            b = self.get(index=index)
            b.name = ""
        else:
            b = BndType(index, "", ())
            self._data.append(b)
        b.name = self._unique_name(name)
        b.color = color
        self._data.sort(key=operator.attrgetter('index'))

    def rem_bnd(self, index):
        for d in self._data:
            if d.index == index:
                self._data.remove(d)
                return

    def get_names(self):
        return [n.name for n in self._data]

    def get(self, name=None, index=None, orderindex=None):
        '->BndType'
        if name is not None:
            for v in self._data:
                if v.name == name:
                    return v
            else:
                raise KeyError('Unknown boundary name')
        if index is not None:
            for v in self._data:
                if v.index == index:
                    return v
            else:
                raise KeyError('Unknown boundary name')
        if orderindex is not None:
            return self._data[orderindex]

        if len(self._data) > 0:
            return self._data[0]
        else:
            #this could be emitted in transitional program states
            #normally it couldn't
            return BndType(0, "None", (0, 0, 0))

    def xml_save(self, xmlnode):
        'save data to xml node'
        for v in self._data:
            nd = ET.SubElement(xmlnode, "BTYPE")
            nd.attrib["index"] = str(v.index)
            ET.SubElement(nd, "NAME").text = str(v.name)

    def xml_load(self, xmlnode):
        'creates BoundaryTypesList from xml node'
        del self._data[:]
        nds = xmlnode.findall("BTYPE")
        for n in nds:
            ind = int(n.attrib["index"].text)
            nm = n.find("NAME").text
            color = (255, 255, 255)
            self.set_bnd(ind, nm, color)


# class BoundaryTypesView(QtGui.QTreeView):
#     "Widget for boundary list representation"
#     def __init__(self, parent=None):
#         super(BoundaryTypesView, self).__init__(parent)
#         self.setHeaderHidden(True)
#         self.setSelectionBehavior(QtGui.QAbstractItemView.SelectItems)
#         self.header().setResizeMode(
#                 QtGui.QHeaderView.ResizeToContents)
#         self.header().setStretchLastSection(False)

#     def viewOptions(self):
#         'overriden'
#         ret = super(BoundaryTypesView, self).viewOptions()
#         ret.decorationAlignment = QtCore.Qt.AlignHCenter \
#                 | QtCore.Qt.AlignCenter
#         ret.decorationPosition = QtGui.QStyleOptionViewItem.Top
#         return ret

#     def mouseDoubleClickEvent(self, event):
#         'overriden'
#         index = self.indexAt(event.pos())
#         if not index.isValid() or index.column() == 3:
#             return
#         else:
#             #edit existing boundary type
#             bnd_types = globvars.actual_data().boundary_types
#             i = index.internalPointer().data(0)
#             bt = bnd_types.get(index=i)
#             dialog = dlgs.EditBoundaryType(bnd_types, bt, self)
#             if dialog.exec_():
#                 i, nm, col = dialog.ret_value()
#                 if bt.index != i and i in bnd_types._ind_set():
#                     txt = "Boundary with index %i already exists. " % i
#                     txt += "Reset it?"
#                     a = QtGui.QMessageBox.question(None, "Confirmation",
#                           txt, QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
#                     if a == QtGui.QMessageBox.No:
#                         return
#                 com = contcom.EditBoundaryType(bt.index, i, nm, col)
#                 globvars.actual_flow().exec_command(com)

#     def _add_dialog(self):
#         globvars.mainWindow._new_bc()

#     def mousePressEvent(self, event):
#         'overriden'
#         index = self.indexAt(event.pos())
#         if event.button() == QtCore.Qt.LeftButton:
#             if index.isValid() and index.column() == 3 and index.row() > 0:
#                 #delete
#                 bnd_types = globvars.actual_data().boundary_types
#                 i = index.internalPointer().data(0)
#                 bt = bnd_types.get(index=i)
#                 com = contcom.EditBoundaryType(bt.index, None, None, None)
#                 globvars.actual_flow().exec_command(com)

#         if event.button() == QtCore.Qt.RightButton:
#             #popup menu in nonactive zone
#             if not index.isValid():
#                 def trig(m, name, fun):
#                     act = QtGui.QAction(name, self)
#                     act.triggered.connect(fun)
#                     m.addAction(act)
#                 menu = QtGui.QMenu(self)
#                 trig(menu, "Add boundary", self._add_dialog)
#                 menu.popup(self.viewport().mapToGlobal(event.pos()))


# class BoundaryTypesModel(qtbp.RowsTreeModel):
#     'table model for boundaries representation'
#     def __init__(self, bt, parent=None):
#         super(BoundaryTypesModel, self).__init__(4, parent)
#         self.bt = bt

#     def reset(self):
#         "overriden"
#         #remove all data from tree representation
#         #and rebuild data
#         self._root_item.clear_childs()
#         self._fill_data_tree(self.bt)
#         super(BoundaryTypesModel, self).reset()

#     def data(self, index, role):
#         "overriden"
#         item = index.internalPointer()
#         if role in [QtCore.Qt.DisplayRole, QtCore.Qt.EditRole]:
#             if index.column() == 2:
#                 return None
#             else:
#                 return item.data(index.column())
#         elif role in [QtCore.Qt.BackgroundColorRole]:
#             if index.column() == 2:
#                 return QtGui.QColor(*item.data(2))
#         elif role in [QtCore.Qt.SizeHintRole]:
#             if index.column() in [2, 3]:
#                 return QtCore.QSize(20, 18)
#         elif role in [QtCore.Qt.DecorationRole]:
#             if index.column() == 3 and index.row() > 0:
#                 return qtbp.get_icon('del')

#     def flags(self, index):
#         "overriden"
#         return QtCore.Qt.ItemIsEnabled

#     def _fill_data_tree(self, bt):
#         names = self.bt.get_names()
#         for n in names:
#             v = self.bt.get(name=n)
#             v = (v.index, v.name, v.color, None)
#             qtbp.TreeItem(v, self._root_item)
