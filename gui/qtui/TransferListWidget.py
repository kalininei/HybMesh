from PyQt4.QtGui import (QWidget, QLineEdit)
import ui_TransferListWidget


class TransferListWidget(QWidget, ui_TransferListWidget.Ui_Form):
    """ Widget with two lists and two buttons
    which are used to transfer list strings from
    one to another
    """

    def __init__(self, parent=None):
        super(TransferListWidget, self).__init__(parent)
        self.setupUi(self)
        self.bt_left.clicked.connect(self._bleft)
        self.bt_right.clicked.connect(self._bright)

    def set_lists(self, left_list, right_list):
        self.set_left_list(left_list)
        self.set_right_list(right_list)

    def set_left_list(self, lst):
        self._set_list_info(self.lst_left, lst)

    def set_right_list(self, lst):
        self._set_list_info(self.lst_right, lst)

    def get_lists(self):
        " -> [list of strings], [list of strings] "
        return self.get_left_list(), self.get_right_list()

    def get_left_list(self):
        " -> [list of strings] "
        return self._get_list_info(self.lst_left)

    def get_right_list(self):
        " -> [list of strings] "
        return self._get_list_info(self.lst_right)

    def _get_list_info(self, lst):
        src = []
        for i in range(lst.count()):
            it = lst.item(i)
            src.append(str(it.text()))
        return src

    def _set_list_info(self, lst, data):
        lst.clear()
        for s in data:
            lst.addItem(s)

    def _bleft(self):
        for it in self.lst_right.selectedItems():
            self.lst_left.addItem(it.text())
            #self.lst_right.removeItemWidget(it)
            self.lst_right.takeItem(self.lst_right.row(it))

    def _bright(self):
        for it in self.lst_left.selectedItems():
            self.lst_right.addItem(it.text())
            #self.lst_left.removeItemWidget(it)
            self.lst_left.takeItem(self.lst_left.row(it))

if __name__ == "__main__":
    import sys
    import PyQt4.QtGui

    app = PyQt4.QtGui.QApplication(sys.argv)
    window = QTransferListWidget()
    window.show()
    sys.exit(app.exec_())



