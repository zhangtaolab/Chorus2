from __future__ import print_function
import sys
import subprocess
import signal
import time
import os
import shlex
from Choruslib.Blatres import Blatres


def start_gfServer(file2bit, gfspath, stepsize=7, blatport=10010):

    '''
    :param file2bit:
    :param gfspath:
    :param stepsize:
    :param blatport:
    :return: blat gfserver pid
    '''

    gfserverpath = '"' + os.path.realpath(gfspath) + '"'

    genomefile = os.path.realpath(file2bit)

    genomefilename = os.path.basename(genomefile)

    genomefilepath = os.path.dirname(genomefile)

    blatcmd = gfserverpath + " start 127.0.0.1 "\
                  +str(blatport) + " -maxDnaHits=20 -stepSize=" \
                  + str(stepsize) + ' ' + genomefilename
    blatcmd = shlex.split(blatcmd)
    print("start gfServer: ", blatcmd)

    p = subprocess.Popen(blatcmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         cwd=genomefilepath)

    print(p.pid)

    return p.pid


def check_gfServer_running():

    p = subprocess.Popen('ps -A', shell=True, stdout=subprocess.PIPE)

    lines = p.stdout.readlines()

    for line in lines:

        if 'gfServer' in line:

            return True

    return False


def stop_gfServer(p=None):

    if p is not None:

        os.kill(p.pid, signal.SIGTERM)

        time.sleep(5)

    else:

        pids = []
        p = subprocess.Popen('ps -A', shell=True, stdout=subprocess.PIPE)

        lines = p.stdout.readlines()

        for line in lines:

            if b'gfServer' in line:

                pids.append(int(line.split()[0]))

        for pid in pids:

            os.kill(pid,signal.SIGTERM)

            time.sleep(10)


def blat_search_sequence(gfcpath, sequence, blatport=10010, minIdentity=75, file2bit='/'):

    gfclientpath = '"'+os.path.realpath(gfcpath)+'"'

    genomepath = os.path.dirname(file2bit)

    # searchcmd = gfclientpath + " -minIdentity="+str(75)+" -nohead 127.0.0.1 "+str(10010)+" " + "." + " /dev/stdin /dev/stdout"

    queryseq = ">%s\n%s" % (sequence, sequence)

    searchcmd = gfclientpath + " -minIdentity="+str(minIdentity)+" -nohead 127.0.0.1 "+str(blatport)+" " + "." + " /dev/stdin /dev/stdout"

    num = 0

    searchcmd = shlex.split(searchcmd)

    # print(searchcmd)

    p = subprocess.Popen(searchcmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, cwd=genomepath)

    p.stdin.write(queryseq.encode('ascii'))

    p.stdin.close()

    for line in p.stdout:

        if line == b"Output is in /dev/stdout\n":

            continue

        num += 1

    return num




def blat_search(blatpath, sequence, minIdentity, samplename, file2bit='/'):

    pslfile = samplename+'.psl'

    blatcmd = blatpath + ' ' +file2bit + ' ' + sequence + ' ' + " -minIdentity="+str(minIdentity) + ' ' + pslfile

    print("blat", blatcmd)

    p = subprocess.call(blatcmd, shell=True)

    return p


def blat_searchpb(gfcpath, sequence, blatport=10010, minIdentity=75, file2bit='/'):

    '''
    :param gfcpath: gfclient path
    :param sequence: seqeunce or probe
    :param blatport: blat port
    :param minIdentity:
    :param file2bit: genome 2bit file for blat
    :return: blatres class
    '''

    gfclientpath = '"'+os.path.realpath(gfcpath)+'"'

    genomepath = os.path.dirname(file2bit)

    queryseq = ">%s\n%s" % (sequence, sequence)

    searchcmd = gfclientpath + " -minIdentity="+str(minIdentity)+" -nohead 127.0.0.1 "+str(blatport)+" " + "." + " /dev/stdin /dev/stdout"

    num = 0

    searchcmd = shlex.split(searchcmd)

    p = subprocess.Popen(searchcmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, cwd=genomepath)

    p.stdin.write(queryseq.encode('ascii'))

    p.stdin.close()

    blatlines = list()

    for line in p.stdout:

        if line == b"Output is in /dev/stdout\n":

            continue

        else:

            blatlines.append(line.decode("utf-8"))

    blatres = Blatres(seq=sequence, blatlines=blatlines)

    return blatres


def build2bit(fato2bitpath, genomefile):
    pass