import numpy as np
from PyQt5 import QtWidgets, QtGui, QtCore
import sys
from ChorusGUI.test_mplOV import Ui_Form
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class FMPT(QtWidgets.QDialog, Ui_Form):

    def __init__(self, parent=None):

        super(FMPT, self).__init__(parent)
        self.setupUi(self)

        self.x = np.arange(0.0, 5.0, 0.01)

        self.y = np.sin(2*np.pi*self.x) + 0.5*np.random.randn(len(self.x))


        # self.fig = Figure()

        # self.ax1 = self.fig.add_subplot(211)

        # self.ax2 = self.fig.add_subplot(212)

        # self.widgetOV.canvas = FigureCanvas(self.fig)

        self.widgetOV.canvas.ax.plot(self.x, self.y, '-')

        # self.widgetOV.canvas.ax2.plot(self.x, self.y, '-')

        self.widgetOV.canvas.draw()


if __name__ == '__main__':


    app = QtWidgets.QApplication(sys.argv)

    tb = FMPT()

    tb.show()


    sys.exit(app.exec_())