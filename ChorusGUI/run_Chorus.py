from ChorusGUI.Chorus_setup import Ui_Form
import sys
from PyQt5 import QtWidgets, QtGui, QtCore
import os
import time


class ChorusSetup(QtWidgets.QDialog, Ui_Form):

    def __init__(self, parent=None):

        super(ChorusSetup, self).__init__(parent)

        self.setupUi(self)

        # self.progressBar.setVisible(False)
        # self.groupBox_pb.setVisible(False)
        # self.groupBox_bt.setVisible(False)
        self.genomefile = ''

        self.inputfile = ''

        self.pjdir = ''

        self.pushButton_setprojectdir.clicked.connect(self.setProject)

        self.pushButton_setgenomefile.clicked.connect(self.setGenomefile)

        self.pushButton_setinputfile.clicked.connect(self.setInputfile)

        self.pushButton_setgfclient.clicked.connect(self.setgfclient)

        self.pushButton_setgfserver.clicked.connect(self.setgfserver)

        n=100

        self.progressBar.setRange(0, n)

        self.myChorus = StartRun(n=n)

        self.myChorus.notifyProgress.connect(self.onProgress)

        self.pushButton_runchorus.clicked.connect(self.getstart)


    def getstart(self):

        if self.genomefile and self.inputfile and self.pjdir:

            reply = QtWidgets.QMessageBox.question(self, "Message", 'Are you sure to run Chorus?', QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)

            if reply == QtWidgets.QMessageBox.Yes:

                self.myChorus.start()

                self.textEdit.clear()

                self.progressBar.setVisible(True)
        else:

            QtWidgets.QMessageBox.warning(self, 'Warning', 'Please check your project file(s)', QtWidgets.QMessageBox.Ok)

    def onProgress(self, text, i):

        self.textEdit.insertPlainText(text)

        self.progressBar.setValue(i)

    def setProject(self):

        pjdir = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose Project Directory")

        if pjdir:

            self.pjdir = pjdir

            self.lineEdit_projectdir.setText(self.pjdir)

    def setGenomefile(self):

        genomefile, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select genome file")

        if genomefile:

            self.genomefile = genomefile

            self.lineEdit_genomefile.setText(self.genomefile)

    def setInputfile(self):

        inputfile, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select input file")

        if inputfile:

            self.inputfile = inputfile

            self.lineEdit_inputfile.setText(self.inputfile)

    def setgfserver(self):

        gfserverpath, _ = QtWidgets.QFileDialog.getOpenFileName(self, "set blat gfserver")

        if gfserverpath:

            self.gfserverpath = gfserverpath

            self.lineEdit_gfserver.setText(gfserverpath)

    def setgfclient(self):

        gfclientpath, _ = QtWidgets.QFileDialog.getOpenFileName(self, "set blat gfclient")

        if gfclientpath:

            self.gfclientpath = gfclientpath

            self.lineEdit_gfclient.setText(gfclientpath)


class StartRun(QtCore.QThread):

    notifyProgress = QtCore.pyqtSignal(str, int)

    def __init__(self, n):

        super(StartRun, self).__init__()

        self.n = n + 1

        print(self.n)

    def run(self):

        for i in range(self.n):

            sometxt = str(i) + "line\n"

            self.notifyProgress.emit(sometxt, i)

            time.sleep(0.1)


if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)

    chsetup = ChorusSetup()

    chsetup.show()

    sys.exit(app.exec_())
