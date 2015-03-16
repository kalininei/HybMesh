#!/usr/bin/env python
'vtk draw area widget'
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


class VTKWidget(QVTKRenderWindowInteractor):
    """ 2d draw area widget.
    """
    #   Should be created with not None parent. Parent should be the
    #   exact QFrame owner of the widget (QFrame, QSplitter etc)

    #   After invocation of window.show() the VtkWidget.Initialize()
    #   should be called.

    #   If window is QDialog then QDialog.exec_ method should be
    #   reimplemented as following:
    #       QDialog.show()
    #       VTKWidget.Initialize()
    #       return super(<class name>, self).exec_()

    def __init__(self, parent):
        super(VTKWidget, self).__init__(parent)
        #renderer
        self.ren = vtk.vtkRenderer()
        self.GetRenderWindow().AddRenderer(self.ren)
        self.ren.SetBackground(0, 0.3, 0.3)
        #2d interactor style
        self.iren().SetInteractorStyle(vtk.vtkInteractorStyleImage())

    def iren(self):
        return self.GetRenderWindow().GetInteractor()

    def _get_actors(self):
        ret = []
        it = self.ren.GetActors().NewIterator()
        while not it.IsDoneWithTraversal():
            ret.append(it.GetCurrentObject())
            it.GoToNextItem()
        return ret

    def set_actor_list(self, actlist):
        ' sets actors for visualization '
        self.ren.RemoveAllViewProps()
        map(self.ren.AddActor, actlist)

    def remove_actor(self, actor):
        'removes vtkActor from RenderWindow'
        self.ren.RemoveActor(actor)

    def add_actor(self, actor):
        'adds vtkActor to RenderWindow'
        self.ren.AddActor(actor)


if __name__ == "__main__":
    import sys
    from PyQt4 import QtGui

    class Dlg(QtGui.QDialog):
        def __init__(self, parent=None):
            super(Dlg, self).__init__(parent)

            self.win1 = VTKWidget(self)
            self.win2 = VTKWidget(self)

            b = QtGui.QPushButton()
            b.clicked.connect(self.accept)

            layout = QtGui.QHBoxLayout(self)
            layout.addWidget(self.win1)
            layout.addWidget(b)
            layout.addWidget(self.win2)

        def exec_(self):
            self.show()
            self.win1.Initialize()
            self.win2.Initialize()
            return super(Dlg, self).exec_()

    app = QtGui.QApplication(sys.argv)
    mainWindow = Dlg()

    if mainWindow.exec_():
        print "CLOSED"

    sys.exit(app.exec_())
