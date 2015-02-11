# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'GridViewOpt.ui'
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
        Dialog.resize(171, 185)
        self.verticalLayout = QtGui.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.cb_cview = QtGui.QCheckBox(Dialog)
        self.cb_cview.setObjectName(_fromUtf8("cb_cview"))
        self.verticalLayout.addWidget(self.cb_cview)
        self.cb_eview = QtGui.QCheckBox(Dialog)
        self.cb_eview.setChecked(True)
        self.cb_eview.setObjectName(_fromUtf8("cb_eview"))
        self.verticalLayout.addWidget(self.cb_eview)
        self.cb_bview = QtGui.QCheckBox(Dialog)
        self.cb_bview.setObjectName(_fromUtf8("cb_bview"))
        self.verticalLayout.addWidget(self.cb_bview)
        self.cb_nind = QtGui.QCheckBox(Dialog)
        self.cb_nind.setObjectName(_fromUtf8("cb_nind"))
        self.verticalLayout.addWidget(self.cb_nind)
        self.cb_cind = QtGui.QCheckBox(Dialog)
        self.cb_cind.setObjectName(_fromUtf8("cb_cind"))
        self.verticalLayout.addWidget(self.cb_cind)
        self.cb_eind = QtGui.QCheckBox(Dialog)
        self.cb_eind.setObjectName(_fromUtf8("cb_eind"))
        self.verticalLayout.addWidget(self.cb_eind)
        self.buttonBox = QtGui.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setCenterButtons(True)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), Dialog.reject)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), Dialog.accept)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Grid View", None))
        self.cb_cview.setText(_translate("Dialog", "Cells", None))
        self.cb_eview.setText(_translate("Dialog", "Edges", None))
        self.cb_bview.setText(_translate("Dialog", "Boundaries", None))
        self.cb_nind.setText(_translate("Dialog", "Node indicies", None))
        self.cb_cind.setText(_translate("Dialog", "Cell indicies", None))
        self.cb_eind.setText(_translate("Dialog", "Edge indicies", None))

