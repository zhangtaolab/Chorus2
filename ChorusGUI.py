__author__ = 'Forrest'
from ChorusGUI.ChorusUI import Ui_MainWindow
import sys
import os
from PyQt5 import QtCore, QtGui, QtWidgets
import time
import platform
from math import log
from pyfasta import Fasta
from multiprocessing import Pool
from Choruslib import bwa
from Choruslib import jellyfish
from Choruslib import prefilter, primer3_filter
import shutil


class ChorusMW(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):

        super(ChorusMW, self).__init__(parent)

        self.setupUi(self)

        self.setprogramedefault()

        """
            pushButton of program page
        """
        #remove BLAT temporarily
        self.comboBox_program_Aln.removeItem(1)

        self.pushButton_program_JFPath.clicked.connect(self.setJellyfish)

        self.comboBox_program_Aln.currentIndexChanged.connect(self.setAln)

        self.pushButton_program_default.clicked.connect(self.setprogramedefault)

        self.pushButton_program_next.clicked.connect(self.gotopageparameter)

        self.pushButton_program_AlnPath.clicked.connect(self.setAln)

        """
            pushButton of parameter page
        """
        self.spinBox_parameter_threads.setMaximum(56)

        self.pushButton_parameter_genomefile.clicked.connect(self.setGenomefile)

        self.pushButton_parameter_inputfile.clicked.connect(self.setInputfile)

        self.pushButton_parameter_samplename.clicked.connect(self.setSample)

        self.pushButton_parameter_previous.clicked.connect(self.gotopageprogram)

        self.pushButton_parameter_next.clicked.connect(self.gotopagerun)

        self.pushButton_parameter_threadsdefault.clicked.connect(self.setThreadsdefault)

        self.pushButton_parameter_pblength.clicked.connect(self.setProbelengthdefault)

        self.pushButton_parameter_homologydefault.clicked.connect(self.setHomologydefault)

        self.pushButton_parameter_dtmdefault.clicked.connect(self.setdTmdefault)

        self.pushButton_parameter_reset.clicked.connect(self.resetparameter)

        self.comboBox_parameter_rprimer.currentIndexChanged.connect(self.setRprimer)

        """
            pushButton of run page
        """
        self.pushButton_run_previous.clicked.connect(self.gotopageparameter)

        """
            Run Chrous
        """

        self.pushButton_run_run.clicked.connect(self.runChours)

        self.pushButton_run_stop.clicked.connect(self.stopChrous)

        # """
        # stderr stdout
        #
        # """
        #
        # self.nowmessage = Message()
        #
        # self.nowmessage.message.connect(self.updataTextBrowser)
        #
        # sys.stderr  = self.nowmessage
        # sys.stdout = self.nowmessage


    def gotopageprogram(self):

        self.toolBox.setCurrentIndex(0)

    def gotopageparameter(self):

        self.toolBox.setCurrentIndex(1)

    def gotopagerun(self):

        self.toolBox.setCurrentIndex(2)

    def setJellyfish(self):

        jellyfishpath, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select jellyfish File", filter='jellyfish')

        if jellyfishpath:

            jfversion = jellyfish.jfversion(jfbin=jellyfishpath)

            QtWidgets.QMessageBox.information(self, "Jelly Fish", "Jelly Fish Version: %s" % jfversion)

            if jfversion == 'None':

                self.lineEdit_program_JFPath.setText(jellyfishpath)

                self.label_program_JFversion.setText(jfversion)

                self.label_program_JFversion.setStyleSheet("color:red")

            else:
                self.lineEdit_program_JFPath.setText(jellyfishpath)

                self.label_program_JFversion.setText(jfversion)

                self.label_program_JFversion.setStyleSheet("color:green")

                self.jellyfishpath = jellyfishpath

    def setprogramedefault(self):

        machine = platform.machine()

        system = platform.system()

        jellyfishpath = 'bin/jellyfish/' + machine + '-' + system + '/bin/' + 'jellyfish'

        bwapath = 'bin/bwa/' + machine + '-' + system + '/bwa'


        # set default jellyfish path
        if os.path.isfile(jellyfishpath):

            self.jfversion = jellyfish.jfversion(jellyfishpath)

            if not self.jfversion == 'None':
                self.label_program_JFversion.setText(self.jfversion)

                self.label_program_JFversion.setStyleSheet("color:green")

                self.lineEdit_program_JFPath.setText(jellyfishpath)

                self.jellyfishpath = jellyfishpath

        else:

            jellyfishpath = shutil.which("jellyfish")

            if jellyfishpath:

                self.jfversion = jellyfish.jfversion(jellyfishpath)

                if not self.jfversion == 'None':
                    self.label_program_JFversion.setText(self.jfversion)

                    self.label_program_JFversion.setStyleSheet("color:green")

                    self.lineEdit_program_JFPath.setText(jellyfishpath)

                    self.jellyfishpath = jellyfishpath

        # set default bwa path
        if os.path.isfile(bwapath):

            self.bwaversion = bwa.bwaversion(bwapath)

            print(bwapath, self.bwaversion)

            if not self.bwaversion == 'None':
                self.label_program_Alnversion.setText(self.bwaversion)

                self.lineEdit_program_AlnPath.setText(bwapath)

                self.alignerpath = bwapath

                self.aligner = 'BWA'

                self.label_program_Alnversion.setStyleSheet("color:green")
        else:

            bwapath = shutil.which("bwa")

            if bwapath:

                self.bwaversion = bwa.bwaversion(bwapath)

                print(bwapath, self.bwaversion)

                if not self.bwaversion == 'None':
                    self.label_program_Alnversion.setText(self.bwaversion)

                    self.lineEdit_program_AlnPath.setText(bwapath)

                    self.alignerpath = bwapath

                    self.aligner = 'BWA'

                    self.label_program_Alnversion.setStyleSheet("color:green")

    def setAln(self):

        pro = self.comboBox_program_Aln.currentText()

        messtxt = "Select " + pro

        # if os.path.exists(alnpath):

        alignerindex = self.comboBox_program_Aln.currentIndex()

        if alignerindex == 0:

            alnpath, _ = QtWidgets.QFileDialog.getOpenFileName(self, messtxt, filter='bwa')

            self.lineEdit_program_AlnPath.setText(alnpath)

            self.bwaversion = bwa.bwaversion(alnpath)

            QtWidgets.QMessageBox.information(self, "BWA", "BWA Version: %s" % self.bwaversion)

            if self.bwaversion == 'None':

                self.label_program_Alnversion.setStyleSheet("color:red")

                self.label_program_Alnversion.setText(self.bwaversion)

            else:

                self.label_program_Alnversion.setText(self.bwaversion)

                self.lineEdit_program_AlnPath.setText(alnpath)

                self.alignerpath = alnpath

                self.aligner = 'BWA'

                self.label_program_Alnversion.setStyleSheet("color:green")

        if alignerindex == 1:

            self.aligner = 'BLAT'

    """
        parameter page
    """

    def setGenomefile(self):

        genomefile, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select Genome File")

        if genomefile:
            self.lineEdit_parameter_genomfile.setText(genomefile)

            # print(jellyfishpath)

    def setInputfile(self):

        inputfile, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select input File")

        if inputfile:
            self.lineEdit_parameter_inputfile.setText(inputfile)

    def setSample(self):

        sampledir = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose a Directory for saving result")

        if sampledir:
            self.lineEdit_parameter_samplefolder.setText(sampledir)

    def setThreadsdefault(self):

        self.spinBox_parameter_threads.setValue(4)

    def setProbelengthdefault(self):

        self.spinBox_parameter_pblength.setValue(45)

    def setHomologydefault(self):

        self.spinBox_parameter_homolog.setValue(75)

    def setdTmdefault(self):

        self.spinBox_parameter_dTm.setValue(10)

    def resetparameter(self):

        self.spinBox_parameter_threads.setValue(4)

        self.spinBox_parameter_pblength.setValue(45)

        self.spinBox_parameter_homolog.setValue(75)

        self.spinBox_parameter_dTm.setValue(10)

    def setRprimer(self):

        currentIndex = self.comboBox_parameter_rprimer.currentIndex()

        if currentIndex == 0:
            self.lineEdit_parameter_rprimer.setText('')
            self.lineEdit_parameter_rprimer.setReadOnly(True)

        if currentIndex == 1:
            self.lineEdit_parameter_rprimer.setText('CGTGGTCGCGTCTCA')
            self.lineEdit_parameter_rprimer.setReadOnly(True)

        if currentIndex == 2:
            self.lineEdit_parameter_rprimer.setReadOnly(False)

    """
        run page
    """


    def runChours(self):

        # self.kmer = int(int(self.spinBox_parameter_pblength.value())/100*int(self.spinBox_parameter_homolog.value())-1)



        self.inputfile = self.lineEdit_parameter_inputfile.text()

        self.genomefile = self.lineEdit_parameter_genomfile.text()

        self.samplefolder = self.lineEdit_parameter_samplefolder.text()

        self.genomesize = int(os.path.getsize(self.lineEdit_parameter_genomfile.text()) / 1000000)

        # use genome size to estimate kmer
        self.kmer = int(log(self.genomesize, 4) + 1)

        if self.kmer < 17:
            self.kmer = 17

        self.jfkmerfile = os.path.join(self.lineEdit_parameter_samplefolder.text(),
                                       (os.path.basename(self.lineEdit_parameter_genomfile.text())
                                        + '_' + str(self.kmer) + 'mer.jf'))

        gsmessage = 'Genome Size: ' + str(self.genomesize) + 'M'

        self.textBrowser_run.append(gsmessage)

        if os.path.isfile(self.jfkmerfile):
            print('find:', self.jfkmerfile)
            kmmess = "Found kmerfile " + self.jfkmerfile + ". Do you want rebuild it?"
            kmbuild = QtWidgets.QMessageBox.question(self, 'Message', kmmess,
                                                     QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                                     QtWidgets.QMessageBox.No)

            if kmbuild == QtWidgets.QMessageBox.No:

                self.kmerbuild = False

            else:

                self.kmerbuild = True

        else:
            print("not find", self.jfkmerfile)
            self.kmerbuild = True

        if self.aligner == 'BWA':

            self.indexfile = os.path.basename(self.genomefile)

            bwatestindex = os.path.join(self.samplefolder, self.indexfile + '.sa')

            if os.path.isfile(bwatestindex):
                print('find:', bwatestindex)
                bwamess = "Found bwa index file " + bwatestindex + ". Do you want rebuild it?"
                bwabuild = QtWidgets.QMessageBox.question(self, 'Message', bwamess,
                                                          QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                                          QtWidgets.QMessageBox.No)

                if bwabuild == QtWidgets.QMessageBox.No:

                    self.indexbuild = False

                else:

                    self.indexbuild = True

            else:
                print("not find", bwatestindex)
                self.indexbuild = True

        self.pool = Pool(self.spinBox_parameter_threads.value())

        self.choursRuner = RunChorus(jellyfishpath=self.jellyfishpath,
                                     alnpath=self.lineEdit_program_AlnPath.text(),
                                     threadsnumber=self.spinBox_parameter_threads.value(),
                                     kmer=self.kmer,
                                     homology=self.spinBox_parameter_homolog.value(),
                                     pblength=self.spinBox_parameter_pblength.value(),
                                     dTm=self.spinBox_parameter_dTm.value(),
                                     inputfile=self.inputfile,
                                     genomefile=self.genomefile,
                                     samplefolder=self.samplefolder,
                                     jfkmerfile=self.jfkmerfile,
                                     indexfile=self.indexfile,
                                     rprimer=self.lineEdit_parameter_rprimer.text(),
                                     aligner=self.aligner,
                                     kmerbuild=self.kmerbuild,
                                     indexbuild=self.indexbuild,
                                     pool=self.pool,
                                     )

        self.choursRuner.notifyProgress.connect(self.updataProgressbar)

        self.choursRuner.notifyMessage.connect(self.updataTextBrowser)

        self.choursRuner.start()

    def updataProgressbar(self, n):

        self.progressBar_run.setValue(int(n))

        self.label_run_percentage.setText(str(n))

    def updataTextBrowser(self, infortext):

        self.textBrowser_run.append(time.strftime('%Y-%m-%d %H:%M:%S'))

        self.textBrowser_run.append(str(infortext))

    def stopChrous(self):

        onstopreply = QtWidgets.QMessageBox.question(self, 'Message', "Are you sure to Stop?",
                                                     QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                                     QtWidgets.QMessageBox.No)
        if onstopreply == QtWidgets.QMessageBox.Yes:

            if self.choursRuner.isRunning():

                self.choursRuner.terminate()

                self.pool.terminate()

                # jellyfish.stop_jellyfish()

                # bwa.stop_bwa()

                self.updataTextBrowser('Process Stopped!!')




class RunChorus(QtCore.QThread):
    notifyProgress = QtCore.pyqtSignal(float)

    notifyMessage = QtCore.pyqtSignal(str)

    finishsingle = QtCore.pyqtSignal(str)

    def __init__(self, jellyfishpath, alnpath, threadsnumber, kmer, homology, pblength, dTm, inputfile, genomefile,
                 samplefolder, jfkmerfile, rprimer, aligner, kmerbuild, indexfile, indexbuild, pool, bfsize='1G', size='100M',
                 lowercount=2):

        super(RunChorus, self).__init__()

        self.jellyfishpath = jellyfishpath

        self.alnpath = os.path.abspath(alnpath)

        self.threadsnumber = threadsnumber

        self.homology = homology

        self.pblength = pblength

        self.dTm = dTm

        self.inputfile = inputfile

        self.genomefile = genomefile

        self.samplefolder = samplefolder

        self.bfsize = bfsize

        self.kmer = kmer

        self.aligner = aligner

        self.indexbuild = indexbuild

        # jfoutfilename = os.path.basename(genomefile) + '_' + str(kmer)+'mer'

        self.jfkmerfile = jfkmerfile

        self.lowercount = lowercount

        self.size = size

        self.rprimer = str(rprimer)

        self.pool = pool

        print(jellyfishpath, alnpath, threadsnumber, kmer, homology, pblength, dTm, inputfile, genomefile,
              samplefolder, jfkmerfile, rprimer, sep='\n')

        # self.pool = QtCore.QThreadPool()
        #
        # self.pool.setMaxThreadCount(self.threadsnumber)

        self.indexfile = indexfile

        self.keepprobe = list()

        self.tmpfilename = os.path.join(samplefolder, 'oligo_tmp.fa')

        self.tmpsam = os.path.join(samplefolder, 'oligo_tmp.sam')

        self.kmerbuild = kmerbuild

        #slide step length for next oligo if pass kmer count
        self.step = 4

        self.progressnumber = 0


    def run(self):

        if self.kmerbuild:

            jfcounter = jellyfish.jfcount(jfpath=self.jellyfishpath, mer=self.kmer,
                                          infile=self.genomefile, output=self.jfkmerfile, threads=self.threadsnumber,
                                          lowercount=self.lowercount, size=self.size)

            """
                check jelly fish count run correctly
            """
            if jfcounter:

                self.progressnumber = self.progressnumber + 5

                self.notifyProgress.emit(self.progressnumber)

                self.notifyMessage.emit("JellyFish Count finished...")

            else:

                self.notifyMessage.emit("JellyFish Count Error!!!")

        else:
            jfcountmess = "Use " + self.jfkmerfile
            self.progressnumber = self.progressnumber + 5
            self.notifyProgress.emit(self.progressnumber)
            self.notifyMessage.emit(jfcountmess)

        if self.indexbuild:

            if self.aligner == 'BWA':

                bwa.bwaindex(self.alnpath, self.genomefile, self.samplefolder)

                self.notifyMessage.emit("BWA Index build finished...")

                self.progressnumber = self.progressnumber + 5
                self.notifyProgress.emit(self.progressnumber)

            elif self.aligner == 'BLAT':

                """
                    add code for BLAT
                """

                pass
        else:

            self.progressnumber = self.progressnumber + 5
            self.notifyProgress.emit(self.progressnumber)

        """
            load and splite input file
        """

        # splite sequence longer than 10M
        spsize = 10000000

        maxkmerscore = int(self.pblength * self.homology / 100) - self.kmer

        jffilteredprobe = list()

        fastain = Fasta(self.inputfile)

        jffpbrunerlist = list()


        for seqname in fastain.keys():

            chrlen = len(fastain[seqname])

            if chrlen < spsize:

                start = 0

                end = chrlen - 1

                jffpbruner = jellyfish.JFfpbruner(jfpath=self.jellyfishpath, jfkmerfile=self.jfkmerfile, mer=self.kmer,
                                                  pyfasta=fastain, seqname=seqname, pblength=self.pblength,
                                                  maxkmerscore=maxkmerscore, start=start,
                                                  end=end, step=self.step)
                jffpbrunerlist.append(jffpbruner)

            else:

                chrblock = int(chrlen / spsize) + 1

                for i in range(chrblock):

                    start = i * spsize

                    end = start + spsize - 1

                    if end >= chrlen:

                        end = chrlen - 1

                    jffpbruner = jellyfish.JFfpbruner(jfpath=self.jellyfishpath, jfkmerfile=self.jfkmerfile, mer=self.kmer,
                                                  pyfasta=fastain, seqname=seqname, pblength=self.pblength,
                                                  maxkmerscore=maxkmerscore, start=start,
                                                  end=end, step=self.step)

                    jffpbrunerlist.append(jffpbruner)



        jffinished = 0

        for curpblist in self.pool.imap_unordered(jellyfish.kmerfilterprobe, jffpbrunerlist):

            jffilteredprobe.extend(curpblist)

            tmpprogress = float(format(self.progressnumber + (jffinished/len(jffpbrunerlist) * 40),".2f"))

            self.notifyProgress.emit(tmpprogress)

            if self.isRunning():

                print("running")

            else:

                print("not running")

            jffinished += 1


        self.notifyMessage.emit('kmer filter finished!!')

        self.progressnumber = 50.0

        self.notifyProgress.emit(self.progressnumber)

        tmppbfa = os.path.join(self.samplefolder, os.path.basename(self.inputfile)+'_tmp_probes.fa')

        tmppbfaio = open(tmppbfa, 'w')

        seqnum = 0

        for tmppb in jffilteredprobe:

            print('>','seq',seqnum, sep='',file=tmppbfaio)


            print(tmppb,file=tmppbfaio)


            seqnum += 1

        tmppbfaio.close()

        #delete jffilteredprobe and release memory
        del jffilteredprobe

        bwaindexfile = os.path.join(self.samplefolder, os.path.basename(self.genomefile))

        bwafiltedpb = bwa.bwafilter(bwabin=self.alnpath, reffile=bwaindexfile, inputfile=tmppbfa, minas=self.pblength,
                                    maxxs=int(self.pblength * self.homology / 100), threadnumber=self.threadsnumber)


        tmpbwaftlist = os.path.join(self.samplefolder, os.path.basename(self.inputfile)+'.bed')

        alltmpbwaftlist = os.path.join(self.samplefolder, os.path.basename(self.inputfile)+'_all.bed')

        tmpbwaftlistio = open(tmpbwaftlist,'w')

        allbwaftlistio = open(alltmpbwaftlist,'w')

        seqlenfile = os.path.join(self.samplefolder, os.path.basename(self.inputfile))+'.len'

        seqlenio = open(seqlenfile, 'w')

        seqlength = bwa.bwareflength(bwabin=self.alnpath, reffile=bwaindexfile)

        for seqname in seqlength:

            print(seqname, seqlength[seqname], sep='\t', file=seqlenio)

        seqlenio.close()


        oligobefortmf = list()

        for pbtmp in bwafiltedpb:

            # print(pbtmp, file=tmpbwaftlistio)
            nowpbcounter = dict()

            nowpbcounter['seq'] = pbtmp

            nowpbcounter['dTm'] = self.dTm

            nowpbcounter['rprimer'] = self.rprimer


            oligobefortmf.append(nowpbcounter)

        keepedprobe = list()

        self.progressnumber = 55

        self.notifyProgress.emit(self.progressnumber)

        ctedpb = 0



        oligobefortmflen = len(oligobefortmf)

        for (pb, keep) in self.pool.imap_unordered(probefilter, oligobefortmf):

            if keep:

                keepedprobe.append(pb)
                # print(pb, file=tmpbwaftlistio)

            ctedpb += 1

            if ctedpb % 10000 == 0:

                tmpprogress = float(format(self.progressnumber + (ctedpb/oligobefortmflen * 30),".2f"))

                self.notifyProgress.emit(tmpprogress)

        self.notifyProgress.emit(90)

        pbdictbychr = dict()

        #load pb to dict
        for pb in keepedprobe:

            # print(pb, file=tmpbwaftlistio)
            seq, chro, start = pb.split('\t')

            start = int(start)

            if chro in pbdictbychr:

                pbdictbychr[chro][start] = seq

            else:

                pbdictbychr[chro] = dict()



                pbdictbychr[chro][start] = seq


        #get lenth of primer
        lenrprimer = len(self.rprimer)

        if lenrprimer == 0:

            lenrprimer = 5

        slidwindow = lenrprimer+self.pblength


        for chro in pbdictbychr:

            startn = 0

            for startnow in sorted(pbdictbychr[chro]):

                endnow = startnow + self.pblength - 1

                print(chro, startnow, endnow, pbdictbychr[chro][startnow],file=allbwaftlistio,sep='\t')

                if startnow > startn+slidwindow:

                    #startn = startnow+slidwindow
                    startn = startnow



                    print(chro, startnow, endnow, pbdictbychr[chro][startnow], file=tmpbwaftlistio, sep='\t')


        tmpbwaftlistio.close()

        allbwaftlistio.close()

        #remove temp fasta file
        # os.remove(tmppbfa)

        self.notifyProgress.emit(100)

        self.notifyMessage.emit('all finished!!')



def probefilter(nowpbcounter):

    sequence, seqname, start = nowpbcounter['seq'].split('\t')

    dTm = nowpbcounter['dTm']

    rPrimer = nowpbcounter['rprimer']

    keep = False

    if prefilter.atcg_filter(sequence, 5):

        # print(self.sequence, " not pass atcg")
        pass

    else:
        if rPrimer:

            if primer3_filter.primer3_filter_withRprimer(sequence, rPrimer, dtm=dTm):

                pass

            else:

                keep = True
        else:

            if primer3_filter.primer3_filter(sequence, dtm=dTm):

                pass

            else:

                keep = True

    return (nowpbcounter['seq'], keep)


# class Message(QtCore.QObject):
#
#     message = QtCore.pyqtSignal(str)
#
#     # def __init__(self, parent=None):
#     #
#     #     super(Message, self).__init__(parent)
#
#     def write(self, text):
#
#         self.message.emit(str(text))
#
#
#     def flush(self):
#
#         pass


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)

    chsetup = ChorusMW()

    chsetup.show()

    sys.exit(app.exec_())
