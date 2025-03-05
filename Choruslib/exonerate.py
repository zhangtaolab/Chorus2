from __future__ import print_function
import sys
import subprocess
import signal
import time
import os
import random
import re
import shlex



def start_exonerate_server(espath, esifile, threadnum=24, exonerateport=10005):
    #espath, path of exonerate-server
    esifilepath = os.path.dirname(esifile)
    os.chdir(esifilepath)
    escmd = espath + ' --input ' + esifile + ' --port ' + str(exonerateport) + ' --maxconnections ' + str(threadnum)

    print(escmd)

    p = subprocess.Popen(escmd, shell=True, stderr=subprocess.PIPE)

    return p


def stop_exonerate_server(p=None):

    if p is not None:

        os.kill(p.pid, signal.SIGTERM)

        time.sleep(5)

    else:

        pids = []

        p = subprocess.Popen('ps -A', shell=True, stdout=subprocess.PIPE)

        lines = p.stdout.readlines()

        for line in lines:

            if b'exonerate-server' in line:

                pids.append(int(line.split()[0]))

        for pid in pids:

            os.kill(pid,signal.SIGTERM)

            time.sleep(10)


def check_exonerate_server():

    p = subprocess.Popen('ps -A', shell=True, stdout=subprocess.PIPE)

    lines = p.stdout.readlines()

    for line in lines:

        if 'exonerate-server' in line:

            return True

    return False


def exonerate_search_sequence(ep, sequence, exonerateport=10005, minIdentity=75):
    #ep, exonerate path
    tmpseqfile = str(sequence) + str(random.randint(1, 10000)) + '.fa'

    # tmpseqfiledir = os.path.dirname(tmpseqfile)

    tmpseqfiletr = os.path.join(os.getcwd(), tmpseqfile)

    # print("file:",tmpseqfiletr)

    tmpio = open(tmpseqfiletr, 'w')

    print('>',sequence, sep='', file=tmpio)

    print(sequence, file=tmpio)

    tmpio.close()

    # epcmd = 'exonerate --refine region --model affine:local -s 100 ' + tmpseqfile + ' localhost:'+str(exonerateport) + ' --showalignment no --showcigar no  --showvulgar no --ryo  \'mapping\\t%ql\\t%tal\\t%pi\\t%ps\\t%s\\n\''
    #https://www.ebi.ac.uk/~guy/exonerate/exonerate.man.html
    epcmd = ep+' --refine region --model affine:local -s 100 ' + tmpseqfiletr + ' localhost:'+str(exonerateport) + ' --showalignment no --showcigar no  --showvulgar no --ryo  \'mapping\\t%ql\t%ei\\n\''

    num = 0

    p = subprocess.Popen(epcmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # p.wait(timeout=15)

    time.sleep(1)

    for line in p.stdout:

        if b'mapping\t' in line:

            #print(line)
            line = line.decode('utf-8')

            line = line.rstrip('\n')

            (mp, qlen, ei) = line.split('\t')

            idt = float(ei)/float(qlen)*100

            # print(idt, minIdentity)

            if idt > minIdentity:

                num +=1

    os.remove(tmpseqfiletr)

    return num


def exonerate_search_sequence2(ep, sequencefile, exonerateport=10005, minIdentity=75):

    epcmd = ep+' --refine region --model affine:local -s 100 ' + sequencefile + ' localhost:'+str(exonerateport) + ' --showalignment no --showcigar no  --showvulgar no --ryo  \'mapping\\t%qi\\t%ql\t%ei\n\''

    print(epcmd)

    num = 0

    p = subprocess.Popen(epcmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # p.wait(timeout=15)

    time.sleep(1)

    for line in p.stdout:

        if b'mapping\t' in line:

            #print(line)
            line = line.decode('utf-8')

            line = line.rstrip('\n')

            (mp, seqname, qlen, ei) = line.split('\t')

            idt = float(ei)/float(qlen)*100

            print(line, idt, minIdentity)

            if idt > minIdentity:

                num +=1
    return num

