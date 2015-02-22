# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ScaleDlg.ui'
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

class Ui_scale_dlg(object):
    def setupUi(self, scale_dlg):
        scale_dlg.setObjectName(_fromUtf8("scale_dlg"))
        scale_dlg.resize(359, 333)
        scale_dlg.setAutoFillBackground(False)
        self.gridLayout_2 = QtGui.QGridLayout(scale_dlg)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.transfer_list = TransferListWidget(scale_dlg)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(1)
        sizePolicy.setHeightForWidth(self.transfer_list.sizePolicy().hasHeightForWidth())
        self.transfer_list.setSizePolicy(sizePolicy)
        self.transfer_list.setObjectName(_fromUtf8("transfer_list"))
        self.gridLayout_2.addWidget(self.transfer_list, 0, 0, 1, 2)
        self.frame = QtGui.QFrame(scale_dlg)
        self.frame.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtGui.QFrame.Raised)
        self.frame.setObjectName(_fromUtf8("frame"))
        self.gridLayout = QtGui.QGridLayout(self.frame)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.cb_preserve = QtGui.QCheckBox(self.frame)
        self.cb_preserve.setChecked(True)
        self.cb_preserve.setObjectName(_fromUtf8("cb_preserve"))
        self.gridLayout.addWidget(self.cb_preserve, 0, 0, 1, 2)
        spacerItem = QtGui.QSpacerItem(94, 19, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 0, 2, 1, 2)
        self.label_2 = QtGui.QLabel(self.frame)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.ed_width = QtGui.QLineEdit(self.frame)
        self.ed_width.setEnabled(True)
        self.ed_width.setObjectName(_fromUtf8("ed_width"))
        self.gridLayout.addWidget(self.ed_width, 1, 1, 1, 2)
        self.cb_units = QtGui.QComboBox(self.frame)
        self.cb_units.setObjectName(_fromUtf8("cb_units"))
        self.cb_units.addItem(_fromUtf8(""))
        self.cb_units.addItem(_fromUtf8(""))
        self.gridLayout.addWidget(self.cb_units, 1, 3, 1, 1)
        self.label = QtGui.QLabel(self.frame)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 2, 0, 1, 1)
        self.ed_height = QtGui.QLineEdit(self.frame)
        self.ed_height.setEnabled(False)
        self.ed_height.setObjectName(_fromUtf8("ed_height"))
        self.gridLayout.addWidget(self.ed_height, 2, 1, 1, 2)
        spacerItem1 = QtGui.QSpacerItem(20, 18, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem1, 2, 3, 1, 1)
        self.gridLayout_2.addWidget(self.frame, 1, 0, 1, 1)
        self.groupBox = QtGui.QGroupBox(scale_dlg)
        self.groupBox.setFlat(False)
        self.groupBox.setCheckable(False)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.verticalLayout = QtGui.QVBoxLayout(self.groupBox)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.rb_tr = QtGui.QRadioButton(self.groupBox)
        self.rb_tr.setObjectName(_fromUtf8("rb_tr"))
        self.verticalLayout.addWidget(self.rb_tr)
        self.rb_cnt = QtGui.QRadioButton(self.groupBox)
        self.rb_cnt.setChecked(True)
        self.rb_cnt.setObjectName(_fromUtf8("rb_cnt"))
        self.verticalLayout.addWidget(self.rb_cnt)
        self.rb_bl = QtGui.QRadioButton(self.groupBox)
        self.rb_bl.setObjectName(_fromUtf8("rb_bl"))
        self.verticalLayout.addWidget(self.rb_bl)
        self.gridLayout_2.addWidget(self.groupBox, 1, 1, 1, 1)
        self.buttonBox = QtGui.QDialogButtonBox(scale_dlg)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.gridLayout_2.addWidget(self.buttonBox, 2, 0, 1, 2)

        self.retranslateUi(scale_dlg)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), scale_dlg.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), scale_dlg.reject)
        QtCore.QObject.connect(self.cb_preserve, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.ed_height.setDisabled)
        QtCore.QMetaObject.connectSlotsByName(scale_dlg)

    def retranslateUi(self, scale_dlg):
        scale_dlg.setWindowTitle(_translate("scale_dlg", "Scaling", None))
        self.cb_preserve.setText(_translate("scale_dlg", "Preserve ratio", None))
        self.label_2.setText(_translate("scale_dlg", "Width", None))
        self.ed_width.setText(_translate("scale_dlg", "100", None))
        self.cb_units.setItemText(0, _translate("scale_dlg", "%", None))
        self.cb_units.setItemText(1, _translate("scale_dlg", "length", None))
        self.label.setText(_translate("scale_dlg", "Height", None))
        self.ed_height.setText(_translate("scale_dlg", "100", None))
        self.groupBox.setTitle(_translate("scale_dlg", "Relative point", None))
        self.rb_tr.setText(_translate("scale_dlg", "Top right", None))
        self.rb_cnt.setText(_translate("scale_dlg", "Center", None))
        self.rb_bl.setText(_translate("scale_dlg", "Bottom left", None))

from TransferListWidget import TransferListWidget
