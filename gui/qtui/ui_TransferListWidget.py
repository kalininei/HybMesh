# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'TransferListWidget.ui'
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

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(390, 241)
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.lst_left = QtGui.QListWidget(Form)
        self.lst_left.setObjectName(_fromUtf8("lst_left"))
        self.gridLayout.addWidget(self.lst_left, 0, 0, 5, 1)
        self.bt_right = QtGui.QPushButton(Form)
        self.bt_right.setText(_fromUtf8(""))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/right_arrow.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.bt_right.setIcon(icon)
        self.bt_right.setObjectName(_fromUtf8("bt_right"))
        self.gridLayout.addWidget(self.bt_right, 3, 1, 1, 1)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem, 4, 1, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem1, 0, 1, 1, 1)
        self.lst_right = QtGui.QListWidget(Form)
        self.lst_right.setObjectName(_fromUtf8("lst_right"))
        self.gridLayout.addWidget(self.lst_right, 0, 2, 5, 1)
        self.bt_left = QtGui.QPushButton(Form)
        self.bt_left.setText(_fromUtf8(""))
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/left_arrow.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.bt_left.setIcon(icon1)
        self.bt_left.setObjectName(_fromUtf8("bt_left"))
        self.gridLayout.addWidget(self.bt_left, 2, 1, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Form", None))

import ComGridRes_rc
