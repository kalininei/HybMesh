' commands history window '
from PyQt4.QtGui import (QWidget,
    QTableWidgetItem, QColor,
    QInputDialog)
from PyQt4 import QtCore
import qtui.ui_HistoryWindow


#############  Commands History Window
class HistoryWindow(QWidget, qtui.ui_HistoryWindow.Ui_HistoryViewer):
    def __init__(self, fl_model):
        super(HistoryWindow, self).__init__()
        self.Flows = fl_model
        self.setupUi(self)
        self._connect_all()
        self.update_info()

    def _disconnect_all(self):
        #insert a comment line
        self.tab_coms.itemChanged.disconnect(self.change_comment)
        self.lst_flows.currentItemChanged.disconnect(self.change_current_flow)

    def _connect_all(self):
        #insert a comment line
        self.tab_coms.itemChanged.connect(self.change_comment)
        self.lst_flows.currentItemChanged.connect(self.change_current_flow)

    def update_info(self):
        self._disconnect_all()
        #place existing flows
        self.lst_flows.clear()
        fnames = self.Flows.get_flow_names()
        for fnm in fnames:
            self.lst_flows.addItem(fnm)
        self.lst_flows.setCurrentRow(fnames.index(
            self.Flows.get_actual_flow_name()))
        #table with commands
        self.set_command_table(self.Flows.get_actual_flow())
        self._connect_all()

    def set_command_table(self, com_flow):
        self.tab_coms.clear()
        #get information about the flow
        tab = com_flow.table()

        #fill rows
        self.tab_coms.setRowCount(len(tab))
        lbs = []
        for i, c in enumerate(tab):
            if c.state == 0 or c.state == 2:
                clr = QColor(190, 190, 190)
            else:
                clr = QColor(0, 0, 0)
            if c.isInitial:
                clr = QColor(0, 0, 255)
            #command
            itm1 = QTableWidgetItem(c.strline)
            itm1.setFlags(QtCore.Qt.NoItemFlags)
            itm1.setTextColor(clr)
            self.tab_coms.setItem(i, 0, itm1)
            #comment
            itm2 = QTableWidgetItem(c.comline)
            itm2.setTextColor(clr)
            self.tab_coms.setItem(i, 1, itm2)
            #row label
            if c.isLast:
                lbs.append('>')
            else:
                lbs.append(str(i))
        self.tab_coms.setVerticalHeaderLabels(lbs)

    #buttons
    @QtCore.pyqtSignature("")
    def on_bt_close_clicked(self):
        self.close()

    @QtCore.pyqtSignature("")
    def on_bt_undo_clicked(self):
        self.Flows.get_actual_flow().undo_prev()

    @QtCore.pyqtSignature("")
    def on_bt_redo_clicked(self):
        self.Flows.get_actual_flow().exec_next()

    @QtCore.pyqtSignature("")
    def on_bt_tostart_clicked(self):
        self.Flows.get_actual_flow().to_zero_state()

    @QtCore.pyqtSignature("")
    def on_bt_checkpoint_clicked(self):
        nm, ok = QInputDialog.getText(self, "New Checkpoint",
                "Enter checkpoint name", text="Flow1")
        if ok:
            self.Flows._checkpoint(str(nm))

    @QtCore.pyqtSignature("")
    def on_bt_delflow_clicked(self):
        print "del flow"

    #change the comment line
    def change_comment(self, item):
        #if comment row
        if (item.column() == 1):
            num = item.row()
            c = self.Flows.get_actual_flow()._commands[num]
            c.set_comment(str(item.text()))

    def change_current_flow(self, prv, cur):
        self.Flows.set_actual_flow(str(prv.text()))
        self.update_info()
