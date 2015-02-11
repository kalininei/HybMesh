# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'UniteGridsDlg.ui'
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

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(267, 371)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/mainwin.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        Dialog.setWindowIcon(icon)
        self.buttonBox = QtGui.QDialogButtonBox(Dialog)
        self.buttonBox.setGeometry(QtCore.QRect(100, 340, 156, 25))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.lst_used = QtGui.QListWidget(Dialog)
        self.lst_used.setGeometry(QtCore.QRect(4, 106, 101, 191))
        self.lst_used.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.lst_used.setObjectName(_fromUtf8("lst_used"))
        self.ed_name = QtGui.QLineEdit(Dialog)
        self.ed_name.setGeometry(QtCore.QRect(69, 13, 191, 22))
        self.ed_name.setObjectName(_fromUtf8("ed_name"))
        self.label = QtGui.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(4, 16, 58, 16))
        self.label.setObjectName(_fromUtf8("label"))
        self.lst_unused = QtGui.QListWidget(Dialog)
        self.lst_unused.setGeometry(QtCore.QRect(152, 106, 111, 191))
        self.lst_unused.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.lst_unused.setObjectName(_fromUtf8("lst_unused"))
        self.bt_left = QtGui.QPushButton(Dialog)
        self.bt_left.setGeometry(QtCore.QRect(114, 162, 31, 25))
        self.bt_left.setText(_fromUtf8(""))
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/left_arrow.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.bt_left.setIcon(icon1)
        self.bt_left.setObjectName(_fromUtf8("bt_left"))
        self.bt_right = QtGui.QPushButton(Dialog)
        self.bt_right.setGeometry(QtCore.QRect(114, 190, 31, 25))
        self.bt_right.setText(_fromUtf8(""))
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/right_arrow.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.bt_right.setIcon(icon2)
        self.bt_right.setObjectName(_fromUtf8("bt_right"))
        self.ed_buf = QtGui.QLineEdit(Dialog)
        self.ed_buf.setGeometry(QtCore.QRect(69, 40, 191, 22))
        self.ed_buf.setObjectName(_fromUtf8("ed_buf"))
        self.label_2 = QtGui.QLabel(Dialog)
        self.label_2.setGeometry(QtCore.QRect(4, 40, 61, 16))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.sld_dens = QtGui.QSlider(Dialog)
        self.sld_dens.setGeometry(QtCore.QRect(70, 70, 151, 23))
        self.sld_dens.setMinimum(0)
        self.sld_dens.setMaximum(10)
        self.sld_dens.setPageStep(1)
        self.sld_dens.setProperty("value", 7)
        self.sld_dens.setOrientation(QtCore.Qt.Horizontal)
        self.sld_dens.setTickPosition(QtGui.QSlider.TicksBothSides)
        self.sld_dens.setObjectName(_fromUtf8("sld_dens"))
        self.label_3 = QtGui.QLabel(Dialog)
        self.label_3.setGeometry(QtCore.QRect(5, 72, 61, 16))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.label_4 = QtGui.QLabel(Dialog)
        self.label_4.setGeometry(QtCore.QRect(230, 73, 61, 16))
        self.label_4.setObjectName(_fromUtf8("label_4"))

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), Dialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), Dialog.reject)
        QtCore.QObject.connect(self.sld_dens, QtCore.SIGNAL(_fromUtf8("valueChanged(int)")), self.label_4.setNum)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Unite Grids", None))
        self.ed_name.setText(_translate("Dialog", "UnitedGrid1", None))
        self.label.setText(_translate("Dialog", "Grid name", None))
        self.ed_buf.setText(_translate("Dialog", "1", None))
        self.label_2.setText(_translate("Dialog", "Buffer zone", None))
        self.label_3.setText(_translate("Dialog", "Density", None))
        self.label_4.setText(_translate("Dialog", "7", None))

import ComGridRes_rc
