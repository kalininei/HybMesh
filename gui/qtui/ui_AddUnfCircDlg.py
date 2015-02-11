# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'AddUnfCircDlg.ui'
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

class Ui_add_unf_circ(object):
    def setupUi(self, add_unf_circ):
        add_unf_circ.setObjectName(_fromUtf8("add_unf_circ"))
        add_unf_circ.resize(292, 198)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/mainwin.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        add_unf_circ.setWindowIcon(icon)
        self.label = QtGui.QLabel(add_unf_circ)
        self.label.setGeometry(QtCore.QRect(4, 5, 58, 16))
        self.label.setObjectName(_fromUtf8("label"))
        self.ed_name = QtGui.QLineEdit(add_unf_circ)
        self.ed_name.setGeometry(QtCore.QRect(66, 5, 131, 22))
        self.ed_name.setObjectName(_fromUtf8("ed_name"))
        self.label_2 = QtGui.QLabel(add_unf_circ)
        self.label_2.setGeometry(QtCore.QRect(4, 59, 40, 16))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.ed_rad = QtGui.QLineEdit(add_unf_circ)
        self.ed_rad.setGeometry(QtCore.QRect(66, 59, 95, 22))
        self.ed_rad.setObjectName(_fromUtf8("ed_rad"))
        self.label_3 = QtGui.QLabel(add_unf_circ)
        self.label_3.setGeometry(QtCore.QRect(4, 86, 71, 16))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.ed_na = QtGui.QLineEdit(add_unf_circ)
        self.ed_na.setGeometry(QtCore.QRect(79, 86, 49, 22))
        self.ed_na.setObjectName(_fromUtf8("ed_na"))
        self.label_4 = QtGui.QLabel(add_unf_circ)
        self.label_4.setGeometry(QtCore.QRect(132, 88, 87, 16))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.ed_nr = QtGui.QLineEdit(add_unf_circ)
        self.ed_nr.setGeometry(QtCore.QRect(223, 86, 65, 22))
        self.ed_nr.setObjectName(_fromUtf8("ed_nr"))
        self.label_5 = QtGui.QLabel(add_unf_circ)
        self.label_5.setGeometry(QtCore.QRect(4, 113, 91, 16))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.ed_coef = QtGui.QLineEdit(add_unf_circ)
        self.ed_coef.setGeometry(QtCore.QRect(109, 113, 51, 22))
        self.ed_coef.setObjectName(_fromUtf8("ed_coef"))
        self.cb_trian = QtGui.QCheckBox(add_unf_circ)
        self.cb_trian.setGeometry(QtCore.QRect(4, 140, 124, 21))
        self.cb_trian.setChecked(True)
        self.cb_trian.setObjectName(_fromUtf8("cb_trian"))
        self.buttonBox = QtGui.QDialogButtonBox(add_unf_circ)
        self.buttonBox.setGeometry(QtCore.QRect(132, 166, 156, 25))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.label_6 = QtGui.QLabel(add_unf_circ)
        self.label_6.setGeometry(QtCore.QRect(4, 32, 16, 16))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.ed_x = QtGui.QLineEdit(add_unf_circ)
        self.ed_x.setGeometry(QtCore.QRect(20, 30, 70, 22))
        self.ed_x.setObjectName(_fromUtf8("ed_x"))
        self.label_7 = QtGui.QLabel(add_unf_circ)
        self.label_7.setGeometry(QtCore.QRect(115, 33, 16, 16))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.ed_y = QtGui.QLineEdit(add_unf_circ)
        self.ed_y.setGeometry(QtCore.QRect(130, 31, 95, 22))
        self.ed_y.setObjectName(_fromUtf8("ed_y"))

        self.retranslateUi(add_unf_circ)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), add_unf_circ.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), add_unf_circ.reject)
        QtCore.QMetaObject.connectSlotsByName(add_unf_circ)

    def retranslateUi(self, add_unf_circ):
        add_unf_circ.setWindowTitle(_translate("add_unf_circ", "Add circle uniform grid block", None))
        self.label.setText(_translate("add_unf_circ", "Grid name", None))
        self.ed_name.setText(_translate("add_unf_circ", "CircleGrid1", None))
        self.label_2.setText(_translate("add_unf_circ", "Radius", None))
        self.ed_rad.setText(_translate("add_unf_circ", "1.0", None))
        self.label_3.setText(_translate("add_unf_circ", "Arch partition", None))
        self.ed_na.setText(_translate("add_unf_circ", "10", None))
        self.label_4.setText(_translate("add_unf_circ", "Radius partition", None))
        self.ed_nr.setText(_translate("add_unf_circ", "10", None))
        self.label_5.setText(_translate("add_unf_circ", "Refinement coef", None))
        self.ed_coef.setText(_translate("add_unf_circ", "1.2", None))
        self.cb_trian.setText(_translate("add_unf_circ", "Triangulate center", None))
        self.label_6.setText(_translate("add_unf_circ", "X", None))
        self.ed_x.setText(_translate("add_unf_circ", "0.0", None))
        self.label_7.setText(_translate("add_unf_circ", "Y", None))
        self.ed_y.setText(_translate("add_unf_circ", "0.0", None))

import ComGridRes_rc
