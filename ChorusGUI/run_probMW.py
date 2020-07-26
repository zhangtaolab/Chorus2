import numpy as np
from PyQt5 import QtWidgets, QtGui, QtCore
import sys
# from ChorusGUI.probMW import Ui_MainWindow
# from ChorusGUI.Probes import Probes
from probMW import Ui_MainWindow
from Probes import Probes
from matplotlib.widgets import SpanSelector
import pandas as pd
import os

class DesMainWD(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):

        super(DesMainWD, self).__init__(parent)

        self.nowchr = 'Chromosome'

        self.setupUi(self)

        self.dockWidget_OV.setVisible(False)

        self.actionLoad_probe.triggered.connect(self.select_bedfile)

        self.comboBox_selectchr.activated[str].connect(self.onActionvated)

        self.pushButton_loadchr.clicked.connect(self.draw_graph)

        self.horizontalSlider_start.valueChanged['int'].connect(self.update_graph)

        self.horizontalSlider_end.valueChanged['int'].connect(self.update_graph)

        self.selectedregionlength = 0

        self.probedir = ''

        # self.spinBox_start.valueChanged['int'].connect(self.update_graph)

        # self.spinBox_end.valueChanged['int'].connect(self.update_graph)
        self.pushButton_addpb.clicked.connect(self.add_probes)

        self.pushButton_delete.clicked.connect(self.del_probes)

        self.pushButton_show.clicked.connect(self.draw_overview)

        self.pushButton_projectdir.clicked.connect(self.setProjetDir)

        self.pushButton_probed.clicked.connect(self.select_bedfile)

        self.pushButton_progenome.clicked.connect(self.setGenomefile)

        self.pushButton_proprobedir.clicked.connect(self.setProbeDir)

        self.pushButton_saveprobe.clicked.connect(self.saveProbe)

        self.sortedperkbcount = object()

    def select_bedfile(self):

        file, _ = QtWidgets.QFileDialog.getOpenFileName()

        print(file)

        if file:

            # self.lineEdit.setText(str(file[0]))

            chrlist = self.get_chr(file)

            self.comboBox_selectchr.addItems(chrlist)

            self.comboBox_selectchr.setCurrentIndex(0)

            self.nowchr = chrlist[0]

            #str(file[0])

            self.label_filename.setText(file)

            self.probeset = self.probe.probe

            self.max_range = int(self.probe.maxlength/1000) + 1

            self.spinBox_end.setMaximum(self.max_range)

            self.horizontalSlider_end.setMaximum(self.max_range)

            self.spinBox_start.setMaximum(self.max_range)

            self.horizontalSlider_start.setMaximum(self.max_range)

            self.label_bedfile.setText(file)

            # self.label.setText(self.nowchr)

    def get_chr(self, filename):

        self.probe = Probes(filename)

        chrlist = self.probe.chrs

        return chrlist

    def onActionvated(self, text):

        self.statusbar.showMessage(text)

        self.nowchr = text

    def draw_graph(self):

        self.nowprobe = self.probeset[self.probeset[0] == self.nowchr]

        self.sortedprobe = self.nowprobe.sort(columns=1)

        self.perkbprobe = self.sortedprobe[3].value_counts(sort=False)

        self.sortedperkbcount = pd.DataFrame(self.perkbprobe).sort_index()

        self.sortedperkbcount = self.sortedperkbcount.reindex(index=range(0, self.max_range), fill_value=0)
        #
        self.spinBox_end.setMaximum(self.probe.chrlens[self.nowchr])

        self.horizontalSlider_end.setMaximum(self.probe.chrlens[self.nowchr])

        self.spinBox_start.setMaximum(self.probe.chrlens[self.nowchr])

        self.horizontalSlider_start.setMaximum(self.probe.chrlens[self.nowchr])

        self.spinBox_start.setValue(0)

        self.spinBox_end.setValue(self.probe.chrlens[self.nowchr])

        self.horizontalSlider_start.setValue(0)

        self.horizontalSlider_end.setValue(self.probe.chrlens[self.nowchr])

        self.widget.canvas.ax1.clear()

        self.widget.canvas.ax2.clear()

        self.widget.canvas.ax1.plot(pd.rolling_mean(self.sortedperkbcount.Kb,100))

        self.widget.canvas.ax1.set_xlim(0, self.probe.chrlens[self.nowchr])

        self.widget.canvas.ax1.set_title(self.nowchr)

        self.widget.canvas.line2, = self.widget.canvas.ax2.plot(self.sortedperkbcount.Kb)
        # self.widget.canvas.ax2.plot(self.sortedperkbcount.Kb)

        self.widget.canvas.ax2.set_xlim(0, self.probe.chrlens[self.nowchr])

        self.widget.canvas.draw()


    def update_graph(self):

        self.widget.canvas.ax2.clear()

        self.widget.canvas.ax2.plot(self.sortedperkbcount.Kb)

        self.widget.canvas.ax2.set_xlim(self.spinBox_start.value(), self.spinBox_end.value())

        self.widget.canvas.ax1.clear()

        self.widget.canvas.ax1.set_title(self.nowchr)

        self.widget.canvas.ax1.plot(pd.rolling_mean(self.sortedperkbcount.Kb,100))

        self.widget.canvas.ax1.axvspan(self.spinBox_start.value(), self.spinBox_end.value(), facecolor=self.comboBox_color.currentText(), alpha=0.5)

        self.widget.canvas.ax1.set_xlim(0, self.probe.chrlens[self.nowchr])

        self.subplotprob = self.nowprobe[self.nowprobe[3] > self.spinBox_start.value()]

        self.subplotprob = self.subplotprob[self.subplotprob[3] < self.spinBox_end.value()]

        self.subplottotalprobe = len(self.subplotprob.index)

        self.horizontalSlider_start.setMaximum(self.spinBox_end.value()-1)

        self.horizontalSlider_end.setMinimum(self.spinBox_start.value()+1)

        self.label_totalpb.setText(str(self.subplottotalprobe))

        self.spinBox_pbnumber.setMaximum(self.subplottotalprobe)

        self.spinBox_pbnumber.setValue(self.subplottotalprobe)

        regionlength = self.horizontalSlider_end.value() - self.horizontalSlider_start.value() + 1

        self.selectedregionlength = regionlength

        mes = "Region Length: "+str(regionlength)+'kb'

        self.statusbar.showMessage(mes)

        self.widget.canvas.draw()

    def oneselect(self, xmins, xmaxs):

        xmins = int(xmins)

        xmaxs = int(xmaxs)

        self.widget.canvas.ax2.clear()

        self.widget.canvas.ax2.plot(self.sortedperkbcount.Kb)

        self.widget.canvas.ax2.set_xlim(xmins, xmaxs)

        self.spinBox_start.setValue(xmins)

        self.spinBox_end.setValue(xmaxs)

        self.subplotprob = self.nowprobe[self.nowprobe[3] < xmaxs]

        self.subplotprob = self.subplotprob[self.subplotprob[3] > xmins]

        self.spinBox_start.setValue(xmins)

        self.spinBox_end.setValue(xmaxs)

        self.subplottotalprobe = len(self.subplotprob.index)

        self.label_totalpb.setText(str(self.subplottotalprobe))

        self.horizontalSlider_start.setMaximum(self.spinBox_end.value()-1)

        self.horizontalSlider_end.setMinimum(self.spinBox_start.value()+1)

        self.spinBox_pbnumber.setMaximum(self.subplottotalprobe)

        self.spinBox_pbnumber.setValue(self.subplottotalprobe)

        # print(self.subplotprob)

        self.widget.canvas.ax1.clear()

        self.widget.canvas.ax1.set_title(self.nowchr)

        self.widget.canvas.ax1.plot(pd.rolling_mean(self.sortedperkbcount.Kb,100))

        self.widget.canvas.ax1.set_xlim(0, self.probe.chrlens[self.nowchr])

        self.widget.canvas.ax1.axvspan(xmins, xmaxs, facecolor=self.comboBox_color.currentText(), alpha=0.5)

        regionlength = self.horizontalSlider_end.value() - self.horizontalSlider_start.value() + 1

        mes = "Region Length: "+str(regionlength)+'kb'

        self.selectedregionlength = regionlength

        self.statusbar.showMessage(mes)

        self.widget.canvas.draw()

    def add_probes(self):

        rowcount = self.tableWidget.rowCount()

        self.tableWidget.insertRow(rowcount)

        #probe density per kb
        pbd = round(self.spinBox_pbnumber.value()/self.selectedregionlength, 1)

        itchr = QtWidgets.QTableWidgetItem(self.nowchr)
        itstart = QtWidgets.QTableWidgetItem(self.spinBox_start.text())
        itend = QtWidgets.QTableWidgetItem(self.spinBox_end.text())
        itcolor = QtWidgets.QTableWidgetItem(self.comboBox_color.currentText())
        ittp = QtWidgets.QTableWidgetItem(self.label_totalpb.text())
        itsp = QtWidgets.QTableWidgetItem(self.spinBox_pbnumber.text())
        itrgl = QtWidgets.QTableWidgetItem(str(self.selectedregionlength))
        itpbd = QtWidgets.QTableWidgetItem(str(pbd))



        qcolor = QtGui.QColor(0,0,0)

        if self.comboBox_color.currentText() == 'green':
            qcolor = QtGui.QColor(0, 255,0)
        if self.comboBox_color.currentText() == 'red':
            qcolor = QtGui.QColor(255, 0,0)

        itcolor.setBackground(qcolor)

        self.tableWidget.setItem(rowcount, 0, itchr)
        self.tableWidget.setItem(rowcount, 1, itstart)
        self.tableWidget.setItem(rowcount, 2, itend)
        self.tableWidget.setItem(rowcount, 3, itcolor)
        self.tableWidget.setItem(rowcount, 4, ittp)
        self.tableWidget.setItem(rowcount, 5, itsp)
        self.tableWidget.setItem(rowcount, 6, itrgl)
        self.tableWidget.setItem(rowcount, 7, itpbd)

    def del_probes(self):

        nowItem = self.tableWidget.currentItem()

        nowit = nowItem.row()

        self.tableWidget.removeRow(nowit)

    def draw_overview(self):

        self.widget_OV.canvas.ax.clear()

        self.widget_OV.canvas.ax.plot(pd.rolling_mean(self.sortedperkbcount.Kb,100))

        rowcount = self.tableWidget.rowCount()

        self.dockWidget_OV.setVisible(True)

        self.widget_OV.canvas.ax.set_title(self.nowchr)

        self.widget_OV.canvas.ax.set_xlim(0, self.probe.chrlens[self.nowchr])

        print("nowchr", self.nowchr)

        for i in range(rowcount):

            itchr = self.tableWidget.item(i, 0).text()

            if itchr == self.nowchr:
                itstart = int(self.tableWidget.item(i,1).text())
                itend = int(self.tableWidget.item(i,2).text())
                itcolor = self.tableWidget.item(i,3).text()

                print(itchr, itstart, itend, itcolor)
                self.widget_OV.canvas.ax.axvspan(itstart, itend, facecolor=itcolor, alpha=0.95)



        regionlength = self.horizontalSlider_end.value() - self.horizontalSlider_start.value() + 1

        self.selectedregionlength = regionlength

        mes = "Region Length: "+str(regionlength)+'kb'

        self.statusbar.showMessage(mes)

        self.widget_OV.canvas.draw()

    def setProjetDir(self):

        # options = QtWidgets.QFileDialog.DontResolveSymlinks | QtWidgets.QFileDialog.ShowDirsOnly

        projectdir = QtWidgets.QFileDialog.getExistingDirectory()

        if projectdir:

            self.projectdir = projectdir

            self.label_prodir.setText(self.projectdir)

    def setGenomefile(self):

        genomefile, _ = QtWidgets.QFileDialog.getOpenFileName()

        if genomefile:

            self.genomefile = genomefile

            self.label_genomefile.setText(self.genomefile)


    def setProbeDir(self):

        probedir = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose Probe Set Directory")

        if probedir:

            self.probedir = probedir

            self.label_probedir.setText(self.probedir)


    def saveProbe(self):

        rowcount = self.tableWidget.rowCount()

        if not self.probedir:

            # self.setProbeDir()

            probedir = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose Probe Set Directory")

            if probedir:

                self.probedir = probedir

                self.label_probedir.setText(self.probedir)

        for i in range(rowcount):

            itchr = self.tableWidget.item(i,0).text()

            itstart = int(self.tableWidget.item(i,1).text())

            itend = int(self.tableWidget.item(i,2).text())

            itcolor = self.tableWidget.item(i,3).text()

            #self.subplotprob = self.nowprobe[self.nowprobe[3] > self.spinBox_start.value()]

            #self.subplotprob = self.subplotprob[self.subplotprob[3] < self.spinBox_end.value()]

            nowprobes = self.probeset[self.probeset[0]==itchr]

            nowprobes = nowprobes[nowprobes[3] > itstart]

            nowprobes = nowprobes[nowprobes[3] < itend]

            nowprobes = nowprobes.drop(3, 1)

            # print(nowprobes)
            outfilename = itcolor + '_' + itchr + '_' + str(itstart) + '_' + str(itend) + '.bed'
            absfile = os.path.join(self.probedir, outfilename)
            nowprobes.to_csv(path_or_buf=absfile, sep='\t', index = False, index_label= False, header=False)




if __name__ == '__main__':

    app = QtWidgets.QApplication(sys.argv)

    tb = DesMainWD()

    tb.show()

    span = SpanSelector(tb.widget.canvas.ax1, tb.oneselect, 'horizontal', useblit=True,
                                 rectprops=dict(alpha=0.3, facecolor='grey'))

    sys.exit(app.exec_())