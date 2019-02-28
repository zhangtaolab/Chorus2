import os
from subprocess import Popen
from subprocess import PIPE
import re
import shutil
from Choruslib import subprocesspath
import time
import signal


def testbwa(bwabin):
    """

    :param bwabin: bwa bin path
    :return: bool, True: bwa tested ok. False: bwa error
    """

    bwacmd = [bwabin]

    bwarun = Popen(bwacmd, stdout=PIPE, stderr=PIPE, shell=True)

    # bwarun.communicate()

    testres = False

    pat = re.compile('Version')

    for i in bwarun.stderr.readlines():

        i = i.decode('utf-8').rstrip('\n')

        if re.search(pat, i):

            testres = True

    bwarun.communicate()

    return testres


def bwaversion(bwabin):
    """

    :param bwabin: bwa bin path
    :return: string, version of bwa
    """

    bwacmd = [bwabin]

    bwarun = Popen(bwacmd, stdout=PIPE, stderr=PIPE)

    pat = re.compile('Version')

    version = 'None'

    for i in bwarun.stderr.readlines():

        i = i.decode('utf-8').rstrip('\n')

        if re.search(pat, i):

            (_, version) = i.split(' ')

    bwarun.communicate()

    return version

def bwaindex(bwabin, reffile, samplefolder):
    """
    bwa index
    :param bwabin: bwa bin path
    :param reffile: reference genome file
    :param samplefolder: sample dir
    :return: no retrun
    """

    refbasename = os.path.basename(os.path.abspath(reffile))

    dscopy = os.path.join(samplefolder, refbasename)

    shutil.copyfile(os.path.abspath(reffile), dscopy)

    # refinsample = os.path.join(samplefolder, refbasename)
    bwabin = os.path.abspath(bwabin)

    bwacmd = [bwabin, 'index', refbasename]
    # print(bwacmd)
    runbwaindex = Popen(bwacmd, cwd=samplefolder)

    runbwaindex.communicate()




def bwaalign(bwabin, reffile, inputfile, outfile, threadnumber=1):
    """
    bwa mem alignment
    :param bwabin: bwa bin path
    :param reffile: reference file, make by bwa index
    :param inputfile: sequence or reads file
    :param outfile: samfile
    :param threadnumber: number of threads
    :return: True
    """

    # bwabin = subprocesspath.subprocesspath(bwabin)

    ##/Users/Forrest/SVN/bwa/bwa mem -O 0 -B 0 -E 0 -k 5 ../DM_404.fa oligo_tmp2.fa
    bwabin = subprocesspath.subprocesspath(bwabin)
    reffile = subprocesspath.subprocesspath(reffile)
    inputfile = subprocesspath.subprocesspath(inputfile)
    outfile = subprocesspath.subprocesspath(outfile)

    bwacmd = ' '.join([bwabin, 'mem', '-O',' 0',' -B',' 0',' -E',' 0',' -k',' 5', '-t',str(threadnumber), reffile, inputfile, '>', outfile])

    print(bwacmd)

    runbwaalign = Popen(bwacmd, shell=True)

    runbwaalign.communicate()

    return True

def samfilter(samfile, minas, maxxs):
    """

    :param samfile: samfile
    :param minas: min AS:i score, suggest probe length
    :param maxxs: max XS:i score, suggest probe length * homology
    :return: list, list of probe/sequence
    """
    seqlist = list()

    pat = re.compile('^@')

    inio = open(samfile,'r')

    aspat = re.compile('AS:i:(\d.)')

    xspat = re.compile('XS:i:(\d.)')

    for i in inio.readlines():

        i = i.rstrip('\n')

        if not re.search(pat, i):

            asmatch = re.search(aspat, i)

            xsmatch = re.search(xspat, i)

            if asmatch:

                asscore = int(asmatch.group(1))

            else:

                continue

            if xsmatch:

                xsscore = int(xsmatch.group(1))

            else:

                continue

            if (asscore >= minas) & (xsscore < maxxs):

                mapinfo = i.split('\t')

                seqlist.append(mapinfo[9])

    return seqlist


def stop_bwa(p=None):

    """
    kill all jellyfish process
    :param p: pid of bwa
    :return: no return
    """

    if p is not None:

        os.kill(p.pid, signal.SIGTERM)

        time.sleep(5)

    else:

        pids = []
        p = Popen('ps -A', shell=True, stdout=PIPE)

        lines = p.stdout.readlines()

        for line in lines:

            if b'bwa' in line:

                pids.append(int(line.split()[0]))

        for pid in pids:

            os.kill(pid,signal.SIGTERM)

            time.sleep(10)


def bwaloci(bwabin, reffile, inputfile, threadnumber=1):

    pat = re.compile('^@')

    bwabin = subprocesspath.subprocesspath(bwabin)
    reffile = subprocesspath.subprocesspath(reffile)
    inputfile = subprocesspath.subprocesspath(inputfile)


    bwacmd = ' '.join([bwabin, 'mem', '-O',' 0',' -B',' 0',' -E',' 0',' -k',' 5', '-t',str(threadnumber), reffile, inputfile])

    print(bwacmd)

    runbwaalign = Popen(bwacmd, shell=True, stdout=PIPE)

    res = list()

    for lin in runbwaalign.stdout.readlines():

        lin = lin.decode('utf-8').rstrip('\n')

        if not re.search(pat, lin):

            infor = lin.split('\t')

            seqnmae = infor[2]

            start = infor[3]

            probeseq = infor[9]

            res.append('\t'.join([probeseq, seqnmae, start]))

    return res


def bwafilter(bwabin, reffile, inputfile, minas, maxxs ,threadnumber=1 ):

    pat = re.compile('^@')

    bwabin = subprocesspath.subprocesspath(bwabin)

    reffile = subprocesspath.subprocesspath(reffile)

    inputfile = subprocesspath.subprocesspath(inputfile)

    bwacmd = ' '.join([bwabin, 'mem', '-O',' 0',' -B',' 0',' -E',' 0',' -k',' 5', '-t',str(threadnumber), reffile, inputfile])

    print(bwacmd)

    aspat = re.compile('AS:i:(\d.)')

    xspat = re.compile('XS:i:(\d.)')

    runbwaalign = Popen(bwacmd, shell=True, stdout=PIPE)

    res = list()

    for lin in runbwaalign.stdout.readlines():
        # print("before decode",lin)
        lin = lin.decode('utf-8').rstrip('\n')
        # print("after decode", lin)
        if not re.search(pat, lin):

            infor = lin.split('\t')

            seqnmae = infor[2]

            start = infor[3]

            probeseq = infor[9]

            asmatch = re.search(aspat, lin)

            xsmatch = re.search(xspat, lin)

            if asmatch:

                asscore = int(asmatch.group(1))

            else:

                continue

            if xsmatch:

                xsscore = int(xsmatch.group(1))

            else:

                continue

            if (asscore >= minas) & (xsscore < maxxs):

                res.append('\t'.join([probeseq, seqnmae, start]))


    runbwaalign.stdout.close()

    runbwaalign.wait()

    return res


    # runbwaalign.communicate()

def bwareflength(bwabin, reffile):

    pat = re.compile('@SQ')

    bwabin = subprocesspath.subprocesspath(bwabin)

    reffile = subprocesspath.subprocesspath(reffile)

    bwacmd = ' '.join([bwabin, 'mem',  reffile, '-'])

    runbwaalign = Popen(bwacmd, shell=True, stdout=PIPE, stdin=PIPE)

    runbwaalign.stdin.write('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'.encode('ascii'))

    runbwaalign.stdin.close()

    seqlength = dict()

    for i in runbwaalign.stdout:

        i = i.decode("utf-8")

        i = i.rstrip('\n')

        if re.search(pat, i):

            (_, seqname, seqlen) = i.split('\t')

            seqname = str(seqname.replace('SN:', ''))

            seqlen = int(seqlen.replace('LN:', ''))

            seqlength[seqname] = seqlen

    return seqlength

    #echo 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' | bwa mem Zea_mays.AGPv3.23.dna.genome.fa -

if __name__ == '__main__':


    bwapath = '../bin/bwa/x86_64-Darwin/bwa'

    seqlength = bwareflength(bwapath, '../Test/DM_404.fa')

    print(seqlength)

    # bwaalign(bwapath, '../Test/DM_404.fa', '../Test/Testsampe/oligo_tmp2.fa', '../Test/Testsampe/outfile.sam',4)
    # bwaindex(bwapath, '../Test/DM_404.fa', '../Test/Testsampe/')

    # bwapath = subprocesspath.subprocesspath(bwapath)
    #

    # seqlist = samfilter('../Test/Testsampe/outfile.sam', minas=45, maxxs=33)
    #
    # for i in seqlist:
    #
    #     print(i)

    # tester = bwaversion(bwapath)
    #
    # print(tester)
    #
    # res = bwaloci(bwapath, '../Test/Testsampe/DM_404.fa', '../Test/Testsampe/DM_test.faprobes.fa',threadnumber=4)
    #
    # for i in res:
    #     print(i)