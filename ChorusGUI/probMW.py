# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'probMW.ui'
#
# Created: Tue Feb 10 10:58:32 2015
#      by: PyQt5 UI code generator 5.4
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1024, 768)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.widget = MplWidget(self.centralwidget)
        self.widget.setObjectName("widget")
        self.verticalLayout.addWidget(self.widget)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.spinBox_start = QtWidgets.QSpinBox(self.centralwidget)
        self.spinBox_start.setMaximum(9999999)
        self.spinBox_start.setObjectName("spinBox_start")
        self.horizontalLayout.addWidget(self.spinBox_start)
        self.horizontalSlider_start = QtWidgets.QSlider(self.centralwidget)
        self.horizontalSlider_start.setMaximum(9999999)
        self.horizontalSlider_start.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_start.setObjectName("horizontalSlider_start")
        self.horizontalLayout.addWidget(self.horizontalSlider_start)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_2.addWidget(self.label_2)
        self.spinBox_end = QtWidgets.QSpinBox(self.centralwidget)
        self.spinBox_end.setMaximum(9999999)
        self.spinBox_end.setProperty("value", 9999999)
        self.spinBox_end.setObjectName("spinBox_end")
        self.horizontalLayout_2.addWidget(self.spinBox_end)
        self.horizontalSlider_end = QtWidgets.QSlider(self.centralwidget)
        self.horizontalSlider_end.setMaximum(9999999)
        self.horizontalSlider_end.setPageStep(10)
        self.horizontalSlider_end.setProperty("value", 9999999)
        self.horizontalSlider_end.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_end.setObjectName("horizontalSlider_end")
        self.horizontalLayout_2.addWidget(self.horizontalSlider_end)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_4.addLayout(self.verticalLayout_2)
        self.verticalLayout_5 = QtWidgets.QVBoxLayout()
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignTop)
        self.groupBox.setFlat(False)
        self.groupBox.setCheckable(False)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.comboBox_color = QtWidgets.QComboBox(self.groupBox)
        self.comboBox_color.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.comboBox_color.setInsertPolicy(QtWidgets.QComboBox.InsertAtBottom)
        self.comboBox_color.setFrame(True)
        self.comboBox_color.setObjectName("comboBox_color")
        self.comboBox_color.addItem("")
        self.comboBox_color.addItem("")
        self.verticalLayout_6.addWidget(self.comboBox_color)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.spinBox_pbnumber = QtWidgets.QSpinBox(self.groupBox)
        self.spinBox_pbnumber.setObjectName("spinBox_pbnumber")
        self.horizontalLayout_6.addWidget(self.spinBox_pbnumber)
        self.label_3 = QtWidgets.QLabel(self.groupBox)
        self.label_3.setObjectName("label_3")
        self.horizontalLayout_6.addWidget(self.label_3)
        self.label_totalpb = QtWidgets.QLabel(self.groupBox)
        self.label_totalpb.setObjectName("label_totalpb")
        self.horizontalLayout_6.addWidget(self.label_totalpb)
        self.verticalLayout_6.addLayout(self.horizontalLayout_6)
        self.pushButton_addpb = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_addpb.setObjectName("pushButton_addpb")
        self.verticalLayout_6.addWidget(self.pushButton_addpb)
        self.verticalLayout_5.addWidget(self.groupBox)
        self.horizontalLayout_4.addLayout(self.verticalLayout_5)
        self.verticalLayout.addLayout(self.horizontalLayout_4)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.label_filename = QtWidgets.QLabel(self.centralwidget)
        self.label_filename.setObjectName("label_filename")
        self.horizontalLayout_5.addWidget(self.label_filename)
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_5.addItem(spacerItem2)
        self.comboBox_selectchr = QtWidgets.QComboBox(self.centralwidget)
        self.comboBox_selectchr.setObjectName("comboBox_selectchr")
        self.horizontalLayout_5.addWidget(self.comboBox_selectchr)
        self.pushButton_loadchr = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_loadchr.setObjectName("pushButton_loadchr")
        self.horizontalLayout_5.addWidget(self.pushButton_loadchr)
        self.verticalLayout.addLayout(self.horizontalLayout_5)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1024, 22))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuShow = QtWidgets.QMenu(self.menubar)
        self.menuShow.setObjectName("menuShow")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.dockWidget_probe = QtWidgets.QDockWidget(MainWindow)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dockWidget_probe.sizePolicy().hasHeightForWidth())
        self.dockWidget_probe.setSizePolicy(sizePolicy)
        self.dockWidget_probe.setObjectName("dockWidget_probe")
        self.dockWidgetContents = QtWidgets.QWidget()
        self.dockWidgetContents.setObjectName("dockWidgetContents")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.dockWidgetContents)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.tableWidget = QtWidgets.QTableWidget(self.dockWidgetContents)
        self.tableWidget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectItems)
        self.tableWidget.setGridStyle(QtCore.Qt.SolidLine)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(8)
        self.tableWidget.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(5, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(6, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(7, item)
        self.tableWidget.horizontalHeader().setDefaultSectionSize(50)
        self.verticalLayout_3.addWidget(self.tableWidget)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem3)
        self.pushButton_saveprobe = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_saveprobe.setObjectName("pushButton_saveprobe")
        self.horizontalLayout_3.addWidget(self.pushButton_saveprobe)
        self.pushButton_show = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_show.setObjectName("pushButton_show")
        self.horizontalLayout_3.addWidget(self.pushButton_show)
        self.pushButton_delete = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_delete.setObjectName("pushButton_delete")
        self.horizontalLayout_3.addWidget(self.pushButton_delete)
        self.verticalLayout_3.addLayout(self.horizontalLayout_3)
        self.verticalLayout_4.addLayout(self.verticalLayout_3)
        self.dockWidget_probe.setWidget(self.dockWidgetContents)
        MainWindow.addDockWidget(QtCore.Qt.DockWidgetArea(2), self.dockWidget_probe)
        self.dockWidget_OV = QtWidgets.QDockWidget(MainWindow)
        self.dockWidget_OV.setEnabled(True)
        self.dockWidget_OV.setFloating(True)
        self.dockWidget_OV.setObjectName("dockWidget_OV")
        self.dockWidgetContents_2 = QtWidgets.QWidget()
        self.dockWidgetContents_2.setObjectName("dockWidgetContents_2")
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout(self.dockWidgetContents_2)
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.widget_OV = MplOVWidget(self.dockWidgetContents_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.widget_OV.sizePolicy().hasHeightForWidth())
        self.widget_OV.setSizePolicy(sizePolicy)
        self.widget_OV.setObjectName("widget_OV")
        self.horizontalLayout_7.addWidget(self.widget_OV)
        self.dockWidget_OV.setWidget(self.dockWidgetContents_2)
        MainWindow.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.dockWidget_OV)
        self.dockWidget_PJ = QtWidgets.QDockWidget(MainWindow)
        self.dockWidget_PJ.setObjectName("dockWidget_PJ")
        self.dockWidgetContents_3 = QtWidgets.QWidget()
        self.dockWidgetContents_3.setObjectName("dockWidgetContents_3")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout(self.dockWidgetContents_3)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.groupBox_2 = QtWidgets.QGroupBox(self.dockWidgetContents_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_2.sizePolicy().hasHeightForWidth())
        self.groupBox_2.setSizePolicy(sizePolicy)
        self.groupBox_2.setTitle("")
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout_8 = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        self.label_4 = QtWidgets.QLabel(self.groupBox_2)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout_8.addWidget(self.label_4)
        self.label_prodir = QtWidgets.QLabel(self.groupBox_2)
        self.label_prodir.setFrameShape(QtWidgets.QFrame.Box)
        self.label_prodir.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.label_prodir.setText("")
        self.label_prodir.setObjectName("label_prodir")
        self.horizontalLayout_8.addWidget(self.label_prodir)
        self.pushButton_projectdir = QtWidgets.QPushButton(self.groupBox_2)
        self.pushButton_projectdir.setObjectName("pushButton_projectdir")
        self.horizontalLayout_8.addWidget(self.pushButton_projectdir)
        self.verticalLayout_8.addLayout(self.horizontalLayout_8)
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.label_5 = QtWidgets.QLabel(self.groupBox_2)
        self.label_5.setObjectName("label_5")
        self.horizontalLayout_9.addWidget(self.label_5)
        self.label_bedfile = QtWidgets.QLabel(self.groupBox_2)
        self.label_bedfile.setMaximumSize(QtCore.QSize(450, 16777215))
        self.label_bedfile.setFrameShape(QtWidgets.QFrame.Box)
        self.label_bedfile.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.label_bedfile.setText("")
        self.label_bedfile.setObjectName("label_bedfile")
        self.horizontalLayout_9.addWidget(self.label_bedfile)
        self.pushButton_probed = QtWidgets.QPushButton(self.groupBox_2)
        self.pushButton_probed.setObjectName("pushButton_probed")
        self.horizontalLayout_9.addWidget(self.pushButton_probed)
        self.verticalLayout_8.addLayout(self.horizontalLayout_9)
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.label_6 = QtWidgets.QLabel(self.groupBox_2)
        self.label_6.setObjectName("label_6")
        self.horizontalLayout_10.addWidget(self.label_6)
        self.label_genomefile = QtWidgets.QLabel(self.groupBox_2)
        self.label_genomefile.setFrameShape(QtWidgets.QFrame.Box)
        self.label_genomefile.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.label_genomefile.setText("")
        self.label_genomefile.setObjectName("label_genomefile")
        self.horizontalLayout_10.addWidget(self.label_genomefile)
        self.pushButton_progenome = QtWidgets.QPushButton(self.groupBox_2)
        self.pushButton_progenome.setObjectName("pushButton_progenome")
        self.horizontalLayout_10.addWidget(self.pushButton_progenome)
        self.verticalLayout_8.addLayout(self.horizontalLayout_10)
        self.horizontalLayout_11 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        self.label_7 = QtWidgets.QLabel(self.groupBox_2)
        self.label_7.setObjectName("label_7")
        self.horizontalLayout_11.addWidget(self.label_7)
        self.label_probedir = QtWidgets.QLabel(self.groupBox_2)
        self.label_probedir.setFrameShape(QtWidgets.QFrame.Box)
        self.label_probedir.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.label_probedir.setText("")
        self.label_probedir.setObjectName("label_probedir")
        self.horizontalLayout_11.addWidget(self.label_probedir)
        self.pushButton_proprobedir = QtWidgets.QPushButton(self.groupBox_2)
        self.pushButton_proprobedir.setObjectName("pushButton_proprobedir")
        self.horizontalLayout_11.addWidget(self.pushButton_proprobedir)
        self.verticalLayout_8.addLayout(self.horizontalLayout_11)
        self.verticalLayout_7.addWidget(self.groupBox_2)
        self.dockWidget_PJ.setWidget(self.dockWidgetContents_3)
        MainWindow.addDockWidget(QtCore.Qt.DockWidgetArea(2), self.dockWidget_PJ)
        self.actionLoad_probe = QtWidgets.QAction(MainWindow)
        self.actionLoad_probe.setObjectName("actionLoad_probe")
        self.actionSave_Project = QtWidgets.QAction(MainWindow)
        self.actionSave_Project.setObjectName("actionSave_Project")
        self.actionQuit = QtWidgets.QAction(MainWindow)
        self.actionQuit.setObjectName("actionQuit")
        self.actionShow_Selected_Probe = QtWidgets.QAction(MainWindow)
        self.actionShow_Selected_Probe.setCheckable(True)
        self.actionShow_Selected_Probe.setObjectName("actionShow_Selected_Probe")
        self.actionExport_Probe = QtWidgets.QAction(MainWindow)
        self.actionExport_Probe.setObjectName("actionExport_Probe")
        self.actionShow_Overview = QtWidgets.QAction(MainWindow)
        self.actionShow_Overview.setCheckable(True)
        self.actionShow_Overview.setChecked(False)
        self.actionShow_Overview.setObjectName("actionShow_Overview")
        self.actionLoad_Project = QtWidgets.QAction(MainWindow)
        self.actionLoad_Project.setObjectName("actionLoad_Project")
        self.actionShow_Project = QtWidgets.QAction(MainWindow)
        self.actionShow_Project.setCheckable(True)
        self.actionShow_Project.setChecked(True)
        self.actionShow_Project.setObjectName("actionShow_Project")
        self.menuFile.addAction(self.actionLoad_probe)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionLoad_Project)
        self.menuFile.addAction(self.actionSave_Project)
        self.menuFile.addAction(self.actionExport_Probe)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionQuit)
        self.menuShow.addAction(self.actionShow_Selected_Probe)
        self.menuShow.addAction(self.actionShow_Overview)
        self.menuShow.addAction(self.actionShow_Project)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuShow.menuAction())

        self.retranslateUi(MainWindow)
        self.horizontalSlider_start.valueChanged['int'].connect(self.spinBox_start.setValue)
        self.spinBox_start.valueChanged['int'].connect(self.horizontalSlider_start.setValue)
        self.spinBox_end.valueChanged['int'].connect(self.horizontalSlider_end.setValue)
        self.horizontalSlider_end.valueChanged['int'].connect(self.spinBox_end.setValue)
        self.dockWidget_probe.visibilityChanged['bool'].connect(self.actionShow_Selected_Probe.setChecked)
        self.actionShow_Selected_Probe.triggered['bool'].connect(self.dockWidget_probe.setVisible)
        self.dockWidget_OV.visibilityChanged['bool'].connect(self.actionShow_Overview.setChecked)
        self.actionShow_Overview.triggered['bool'].connect(self.dockWidget_OV.setVisible)
        self.dockWidget_PJ.visibilityChanged['bool'].connect(self.actionShow_Project.setChecked)
        self.actionShow_Project.triggered['bool'].connect(self.dockWidget_PJ.setVisible)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Chrous"))
        self.label.setText(_translate("MainWindow", "Start"))
        self.label_2.setText(_translate("MainWindow", "End "))
        self.groupBox.setTitle(_translate("MainWindow", "Choose Probe Color"))
        self.comboBox_color.setItemText(0, _translate("MainWindow", "red"))
        self.comboBox_color.setItemText(1, _translate("MainWindow", "green"))
        self.label_3.setText(_translate("MainWindow", "/"))
        self.label_totalpb.setText(_translate("MainWindow", "Totalpb"))
        self.pushButton_addpb.setText(_translate("MainWindow", "Add"))
        self.label_filename.setText(_translate("MainWindow", "file name"))
        self.pushButton_loadchr.setText(_translate("MainWindow", "Load Chr"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuShow.setTitle(_translate("MainWindow", "Show"))
        self.dockWidget_probe.setWindowTitle(_translate("MainWindow", "Probes"))
        self.tableWidget.setSortingEnabled(True)
        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Chr"))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Start"))
        item = self.tableWidget.horizontalHeaderItem(2)
        item.setText(_translate("MainWindow", "End"))
        item = self.tableWidget.horizontalHeaderItem(3)
        item.setText(_translate("MainWindow", "Color"))
        item = self.tableWidget.horizontalHeaderItem(4)
        item.setText(_translate("MainWindow", "Total Probe"))
        item = self.tableWidget.horizontalHeaderItem(5)
        item.setText(_translate("MainWindow", "Selected Probe"))
        item = self.tableWidget.horizontalHeaderItem(6)
        item.setText(_translate("MainWindow", "Region Length"))
        item = self.tableWidget.horizontalHeaderItem(7)
        item.setText(_translate("MainWindow", "Probe Density"))
        self.pushButton_saveprobe.setText(_translate("MainWindow", "Save"))
        self.pushButton_show.setText(_translate("MainWindow", "Show"))
        self.pushButton_delete.setText(_translate("MainWindow", "Delete"))
        self.dockWidget_OV.setWindowTitle(_translate("MainWindow", "Overview"))
        self.dockWidget_PJ.setWindowTitle(_translate("MainWindow", "Project"))
        self.label_4.setText(_translate("MainWindow", "Project Dir"))
        self.pushButton_projectdir.setText(_translate("MainWindow", "Set"))
        self.label_5.setText(_translate("MainWindow", "BedFile"))
        self.pushButton_probed.setText(_translate("MainWindow", "Set"))
        self.label_6.setText(_translate("MainWindow", "GenomeFile"))
        self.pushButton_progenome.setText(_translate("MainWindow", "Set"))
        self.label_7.setText(_translate("MainWindow", "ProbeDir"))
        self.pushButton_proprobedir.setText(_translate("MainWindow", "Set"))
        self.actionLoad_probe.setText(_translate("MainWindow", "Load Probe"))
        self.actionSave_Project.setText(_translate("MainWindow", "Save Project"))
        self.actionQuit.setText(_translate("MainWindow", "Quit"))
        self.actionShow_Selected_Probe.setText(_translate("MainWindow", "Show Selected Probe"))
        self.actionExport_Probe.setText(_translate("MainWindow", "Export Probe"))
        self.actionShow_Overview.setText(_translate("MainWindow", "Show Overview"))
        self.actionLoad_Project.setText(_translate("MainWindow", "Load Project"))
        self.actionShow_Project.setText(_translate("MainWindow", "Show Project"))

from ChorusGUI.mplwidget import MplWidget
from ChorusGUI.mplOVwidget import MplOVWidget
