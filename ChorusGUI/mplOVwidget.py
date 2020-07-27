import matplotlib
matplotlib.use("Qt5Agg")
from PyQt5 import QtWidgets
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

class MplOVCanvas(FigureCanvas ):

    def __init__(self, subplotnum = 1):
        # self.x = x
        self.fig = Figure()

        # rownumber = int(subplotnum+1/2)


        self.ax = self.fig.add_subplot(111)


        FigureCanvas.__init__(self, self.fig)

        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)

        FigureCanvas.updateGeometry(self)


class MplOVWidget(QtWidgets.QWidget):

    def __init__(self, parent=None):

        QtWidgets.QWidget.__init__(self, parent)

        # self.main_widget = QtWidgets.QWidget(self)
        self.main_fram = QtWidgets.QWidget()

        # self.canvas = MplOVCanvas()
        self.canvas = MplOVCanvas()

        self.canvas.setParent(self.main_fram)

        self.ntb = NavigationToolbar(self.canvas, self.main_fram)

        self.vbl = QtWidgets.QVBoxLayout()

        self.vbl.addWidget(self.ntb)

        self.vbl.addWidget(self.canvas)

        self.setLayout(self.vbl)