import matplotlib
matplotlib.use("Qt5Agg")
from PyQt5 import QtWidgets
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


class MplCanvas(FigureCanvas):

    def __init__(self):
        # self.x = x
        self.fig = Figure()

        self.fig.set_tight_layout(True)

        self.ax1 = self.fig.add_subplot(211)

        self.ax1.set_title("Whole Chr (Mb)")

        self.ax2 = self.fig.add_subplot(212)

        self.ax2.set_title("Zoom in Region (Kb)")

        self.line2 = object()

        FigureCanvas.__init__(self, self.fig)

        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)

        FigureCanvas.updateGeometry(self)


class MplWidget(QtWidgets.QWidget):

    def __init__(self, parent=None):

        QtWidgets.QWidget.__init__(self, parent)

        self.canvas = MplCanvas()

        self.vbl = QtWidgets.QVBoxLayout()

        self.vbl.addWidget(self.canvas)

        self.setLayout(self.vbl)

