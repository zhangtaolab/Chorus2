# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'test_mplOV.ui'
#
# Created: Thu Feb  5 10:35:34 2015
#      by: PyQt5 UI code generator 5.4
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(400, 300)
        self.verticalLayout = QtWidgets.QVBoxLayout(Form)
        self.verticalLayout.setObjectName("verticalLayout")
        self.widgetOV = MplOVWidget(Form)
        self.widgetOV.setObjectName("widgetOV")
        self.verticalLayout.addWidget(self.widgetOV)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "OverView"))

from mplOVwidget import MplOVWidget
