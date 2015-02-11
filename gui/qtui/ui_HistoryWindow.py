# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'HistoryWindow.ui'
#
# Created by: PyQt4 UI code generator 4.11.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_HistoryViewer(object):
    def setupUi(self, HistoryViewer):
        HistoryViewer.setObjectName(_fromUtf8("HistoryViewer"))
        HistoryViewer.resize(638, 378)
        self.bt_close = QtGui.QPushButton(HistoryViewer)
        self.bt_close.setGeometry(QtCore.QRect(540, 340, 90, 31))
        self.bt_close.setObjectName(_fromUtf8("bt_close"))
        self.lst_flows = QtGui.QListWidget(HistoryViewer)
        self.lst_flows.setGeometry(QtCore.QRect(0, 10, 171, 321))
        self.lst_flows.setEditTriggers(QtGui.QAbstractItemView.DoubleClicked|QtGui.QAbstractItemView.EditKeyPressed)
        self.lst_flows.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.lst_flows.setSelectionRectVisible(True)
        self.lst_flows.setObjectName(_fromUtf8("lst_flows"))
        self.tab_coms = QtGui.QTableWidget(HistoryViewer)
        self.tab_coms.setGeometry(QtCore.QRect(180, 10, 451, 321))
        self.tab_coms.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.tab_coms.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.tab_coms.setColumnCount(2)
        self.tab_coms.setObjectName(_fromUtf8("tab_coms"))
        self.tab_coms.setRowCount(0)
        self.bt_undo = QtGui.QPushButton(HistoryViewer)
        self.bt_undo.setGeometry(QtCore.QRect(290, 340, 51, 31))
        self.bt_undo.setObjectName(_fromUtf8("bt_undo"))
        self.bt_redo = QtGui.QPushButton(HistoryViewer)
        self.bt_redo.setGeometry(QtCore.QRect(350, 340, 51, 31))
        self.bt_redo.setObjectName(_fromUtf8("bt_redo"))
        self.bt_tostart = QtGui.QPushButton(HistoryViewer)
        self.bt_tostart.setGeometry(QtCore.QRect(180, 340, 61, 31))
        self.bt_tostart.setObjectName(_fromUtf8("bt_tostart"))
        self.bt_checkpoint = QtGui.QPushButton(HistoryViewer)
        self.bt_checkpoint.setGeometry(QtCore.QRect(0, 342, 81, 31))
        self.bt_checkpoint.setObjectName(_fromUtf8("bt_checkpoint"))
        self.bt_delflow = QtGui.QPushButton(HistoryViewer)
        self.bt_delflow.setGeometry(QtCore.QRect(90, 342, 81, 31))
        self.bt_delflow.setObjectName(_fromUtf8("bt_delflow"))

        self.retranslateUi(HistoryViewer)
        QtCore.QMetaObject.connectSlotsByName(HistoryViewer)

    def retranslateUi(self, HistoryViewer):
        HistoryViewer.setWindowTitle(_translate("HistoryViewer", "Commands history", None))
        self.bt_close.setText(_translate("HistoryViewer", "Close", None))
        self.bt_undo.setText(_translate("HistoryViewer", "undo", None))
        self.bt_redo.setText(_translate("HistoryViewer", "redo", None))
        self.bt_tostart.setText(_translate("HistoryViewer", "to start", None))
        self.bt_checkpoint.setText(_translate("HistoryViewer", "CheckPoint", None))
        self.bt_delflow.setText(_translate("HistoryViewer", "Delete", None))

