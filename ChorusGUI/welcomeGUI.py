# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'welcomeGUI.ui'
#
# Created: Wed Feb  4 09:52:19 2015
#      by: PyQt5 UI code generator 5.4
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(559, 386)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox = QtWidgets.QGroupBox(Dialog)
        self.groupBox.setTitle("")
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.radioButton_new = QtWidgets.QRadioButton(self.groupBox)
        self.radioButton_new.setObjectName("radioButton_new")
        self.verticalLayout_3.addWidget(self.radioButton_new)
        self.radioButton_select = QtWidgets.QRadioButton(self.groupBox)
        self.radioButton_select.setObjectName("radioButton_select")
        self.verticalLayout_3.addWidget(self.radioButton_select)
        self.radioButton_review = QtWidgets.QRadioButton(self.groupBox)
        self.radioButton_review.setObjectName("radioButton_review")
        self.verticalLayout_3.addWidget(self.radioButton_review)
        self.verticalLayout.addWidget(self.groupBox)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.pushButton_OK = QtWidgets.QPushButton(Dialog)
        self.pushButton_OK.setObjectName("pushButton_OK")
        self.horizontalLayout.addWidget(self.pushButton_OK)
        self.verticalLayout_2.addLayout(self.horizontalLayout)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.radioButton_new.setText(_translate("Dialog", "New Project"))
        self.radioButton_select.setText(_translate("Dialog", "Probe Select from Project"))
        self.radioButton_review.setText(_translate("Dialog", "Probe Review from Project"))
        self.pushButton_OK.setText(_translate("Dialog", "OK"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

