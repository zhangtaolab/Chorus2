from os import path
import subprocess
import platform
import os
import signal
import time
from Choruslib import subprocesspath


def sysinfo():
    """

    :return: machine (x86, x86_64, etc) and system (Linux, Darwin, etc) type
    """

    machine = platform.machine()

    system = platform.system()

    return(machine, system)


def jfcount(jfpath, mer, output, infile,threads=1,  size='100M', lowercount=2):
    """
    Only keep >=2 kerm, if kmer==1 score =0
    :param jfpath:
    :param mer:
    :param output:
    :param infile:
    :param threads:
    :param size:
    :param lowercount:
    :return:
    """

    jfpath = subprocesspath.subprocesspath(jfpath)

    output = subprocesspath.subprocesspath(output)

    infile = subprocesspath.subprocesspath(infile)

    jfcountcommand = ' '.join([jfpath, 'count', '--canonical', '-m', str(mer), '-L', str(lowercount),
                               '-t', str(threads), '-o', str(output),  '-s', str(size),  infile])

    print(jfcountcommand)

    p = subprocess.Popen(jfcountcommand, shell=True)

    try:

        outs, errs = p.communicate()

        return True

    except Exception:

        p.kill()

        outs, errs = p.communicate()

        print("Something wrong in jellyfish count")

        return False

#../jellyfish/x86_64-Darwin/bin/jellyfish count  -m 32 -s 100M -t 4 --bc Zea_mays_32mer.bc -L 2 -o Zea_mays_32mer_L2.jf  ../TestData/Zea_mays.AGPv3.23.dna.genome.fa


def jfgeneratorscount(jfpath, mer, output, generators,threads=1,  size='100M'):
    """
    :param jfpath:
    :param mer:
    :param output:
    :param infile:
    :param threads:
    :param size:
    :param lowercount:
    :return:
    """

    jfpath = subprocesspath.subprocesspath(jfpath)

    output = subprocesspath.subprocesspath(output)

    generators = subprocesspath.subprocesspath(generators)

    jfcountcommand = ' '.join([jfpath, 'count', '--canonical', '-m', str(mer), '-g', generators,
                               '-t', str(threads), '-o', str(output),  '-s', str(size)])

    print(jfcountcommand)

    p = subprocess.Popen(jfcountcommand, shell=True)

    try:

        outs, errs = p.communicate()

        return True

    except Exception:

        p.kill()

        outs, errs = p.communicate()

        print("Something wrong in jellyfish count")

        return False


def jfcountbf(jfpath, mer,  output, infile, threads=1,  size='100M', bfsize="1G", lowercount=1):

    '''
    :param jfpath: jellyfish path
    :param mer: kmer
    :param output: output filename
    :param infile: input filename
    :param threads: threads number
    :param size: cache size for jellyfish
    :param bfsize: buffer size for jelly fish
    :param lowercount: filter kmer less than, 2 in default means filter unique kmer
    :return: boolean
    '''

    jfpath = jfpath

    jfcountcommand = ' '.join([jfpath, 'count', '--canonical', '--bf-size', str(bfsize), '-L', str(lowercount),
                               '-m', str(mer), '-t', str(threads), '-o', str(output),  '-s', str(size), infile])

    p = subprocess.Popen(jfcountcommand, shell=True)

    try:

        outs, errs = p.communicate()

        return True

    except Exception:

        p.kill()

        outs, errs = p.communicate()

        print("Something wrong in jellyfish count")

        return False

def jfquery(jfpath, countfile, sequence):
    """

    :param jfpath: jelly fish path
    :param countfile: jellyfish kmer count file
    :param sequence: sequence/kmer sequence
    :return: sequence's kmer count
    """

    # if path.exists(countfile):

    jfquerycommand = ' '.join([jfpath, 'query', countfile, sequence])

    kmerct = subprocess.Popen(jfquerycommand, shell=True, stdout=subprocess.PIPE)

    for lin in kmerct.stdout:

            lin = lin.decode('utf-8')

            lin = lin.rstrip('\n')

            (kmer, number) = lin.split(' ')

            return int(number)

    # else:
    #
    #     print('can\'t open ', countfile)
    #
    #     exit()


def tjf(jfbin):
    """

    :param jfbin: jelly fish bin path
    :return: bool, True: jellyfish works, False: jellyfish error
    """

    tjf = False

    testcmd = jfbin + '  --version'

    res = subprocess.Popen(testcmd, shell=True, stdout=subprocess.PIPE)

    for lin in res.stdout:

        lin = lin.decode('utf-8')

        lin = lin.rstrip('\n')

        (jf, version) = lin.split(' ')

        print(version)

        if version:

            tjf = True

    return tjf


def localjf():

    (machine, system) = sysinfo()

    jfplantform = machine+'-'+system

    jfbin = 'jellyfish/'+jfplantform+ '/bin/'+'jellyfish'

    if not path.isfile(jfbin):

        print("Can not find ", jfbin)

        exit()

    if not tjf(jfbin):

        print("Can not run jellyfish ")

        exit()

    return jfbin


def jfversion(jfbin):

    """

    :param jfbin: jellyfish bin path
    :return: verison of jellyfish, error return 'None'
    """

    # jfbin = path.abspath(jfbin)

    # testcmd = ' '.join([jfbin , '--version'])
    testcmd = [jfbin , '--version']
    # res = subprocess.Popen(testcmd, shell=True, stdout=subprocess.PIPE)
    res = subprocess.Popen(testcmd, stdout=subprocess.PIPE)
    for lin in res.stdout:

        lin = lin.decode('utf-8')

        lin = lin.rstrip('\n')

        (jf, version) = lin.split(' ')

        print(version)

        if version:

            return version

        else:

            return "None"


def jfseqkmerscore(jfpath, jfkmerfile, mer, sequence, maxkmerscore, bfcount=False):
    """

    :param jfpath: use "better path"
    :param jfkmerfile: use "better path"
    :param mer:
    :param sequence:
    :param bfcount:
    :return:
    """

    kmerscore = 0

    if bfcount:

        seqlen = len(sequence)

        mer = int(mer)

        end = mer

        while (end > seqlen):

            start = end - mer

            subseq = sequence[start:end]

            kmerscorenow = jfquery(jfpath=jfpath, countfile=jfkmerfile, sequence=subseq)

            if kmerscorenow > 0:

                kmerscore = kmerscorenow + kmerscore

    else:

        seqlen = len(sequence)

        mer = int(mer)

        end = mer

        while (end <= seqlen):

            start = end - mer

            subseq = sequence[start:end]

            kmerscorenow = jfquery(jfpath=jfpath, countfile=jfkmerfile, sequence=subseq)

            if kmerscorenow > 0:

                # print(sequence, subseq, kmerscorenow)

                kmerscore = kmerscorenow - 1 + kmerscore

            if kmerscore > maxkmerscore:

                break

            # else:
            #
            #     print(sequence, subseq, kmerscorenow)

            end += 1

    return kmerscore


def stop_jellyfish(p=None):

    """
    kill all jellyfish process
    :param p: pid of jellyfish
    :return: no return
    """

    if p is not None:

        os.kill(p.pid, signal.SIGTERM)

        time.sleep(1)

    else:



        pids = []

        p = subprocess.Popen('ps -A', shell=True, stdout=subprocess.PIPE)

        lines = p.stdout.readlines()

        for line in lines:

            if b'jellyfish' in line:

                pids.append(int(line.split()[0]))

        for pid in pids:

            os.kill(pid,signal.SIGTERM)

            # time.sleep(1)



def jfseqkmercount(jfpath, jfkmerfile, mer, sequence, bfcount=False):

    """
    :param jfpath: jellyfish bin path
    :param jfkmerfile: jellyfish kmer count file
    :param mer: int, kmer
    :param sequence: string, sequence for kmerscore count
    :param bfcount:
    :return: list, kmerscore list
    """

    seqlen = len(sequence)

    jfpath = subprocesspath.subprocesspath(jfpath)

    jfkmerfile = subprocesspath.subprocesspath(jfkmerfile)

    jfquerycommand = ' '.join([jfpath, 'query', '-i', '-l', jfkmerfile])

    print(jfquerycommand)

    kmerct = subprocess.Popen(jfquerycommand, shell=True, stdout=subprocess.PIPE,
                              stdin=subprocess.PIPE)


    mer = int(mer)

    end = mer

    jfkmercount = list()

    while (end <= seqlen):

        start = end - mer

        subseq = sequence[start:end]+'\n'

        kmerct.stdin.write(subseq.encode('ascii'))

        kmerct.stdin.flush()

        lin = kmerct.stdout.readline().decode('utf-8').rstrip('\n')

        number = int(lin)

        if number == 2:

            number = 1

        if number > 2:

            number = 2

        jfkmercount.append(number)
        end += 1

    kmerct.stdin.close()

    kmerct.stdout.close()

    # kmerct.terminate()

    kmerct.wait()

    return jfkmercount


def makegenerator(filenames, type='gz',generators='generators'):
    """

    :param filenames:
    :return:
    """

    generatorsio = open(generators,'w')

    if type == 'gz':

        for filename in filenames:

            fileabpath = path.abspath(filename)

            print('gunzip -c', fileabpath, file=generatorsio)

    if type =='text':

        for filename in filenames:

            fileabpath = path.abspath(filename)

            print('cat', fileabpath, file=generatorsio)

    generatorsio.close()


def jfseqkmercountforfilter(jfpath, jfkmerfile, mer, sequence, bfcount=False):

    """
    :param jfpath: jellyfish bin path
    :param jfkmerfile: jellyfish kmer count file
    :param mer: int, kmer
    :param sequence: string, sequence for kmerscore count
    :param bfcount:
    :return: list, kmerscore list
    """

    seqlen = len(sequence)

    jfpath = subprocesspath.subprocesspath(jfpath)

    jfkmerfile = subprocesspath.subprocesspath(jfkmerfile)

    jfquerycommand = ' '.join([jfpath, 'query', '-i', '-l', jfkmerfile])

    print(jfquerycommand)

    kmerct = subprocess.Popen(jfquerycommand, shell=True, stdout=subprocess.PIPE,
                              stdin=subprocess.PIPE)


    mer = int(mer)

    end = mer

    jfkmercount = list()

    while (end <= seqlen):

        start = end - mer

        subseq = sequence[start:end]+'\n'

        kmerct.stdin.write(subseq.encode('ascii'))

        kmerct.stdin.flush()

        lin = kmerct.stdout.readline().decode('utf-8').rstrip('\n')

        number = int(lin)

        jfkmercount.append(number)
        end += 1

    kmerct.stdin.close()

    kmerct.stdout.close()

    # kmerct.terminate()

    kmerct.wait()

    return jfkmercount

# def kmerfilterprobe(jfpath, jfkmerfile, mer, sequence, pblength, maxkmerscore, step=4, bfcount=False):
def kmerfilterprobe(jffpbruner):
    """

    :param jffpbruner: JFfpbrunner class
    :return: probe list
    """
    sequence = jffpbruner.pyfasta[jffpbruner.seqname][jffpbruner.start:jffpbruner.end]

    print("kmer filtering ", jffpbruner.seqname+":"+str(jffpbruner.start)+':'+str(jffpbruner.end))

    jfkmercount = jfseqkmercount(jffpbruner.jfpath, jffpbruner.jfkmerfile, jffpbruner.mer, sequence)

    st = 0

    seqlength = len(sequence)

    probelist = list()

    while (st < (seqlength-jffpbruner.pblength)):

        sp = st + jffpbruner.pblength

        nowpb = sequence[st:sp]

        kmerscore = sum(jfkmercount[st:(sp - 1)])

        st += 1

        if kmerscore < jffpbruner.maxkmerscore:

            probelist.append(nowpb)

            st += jffpbruner.step

    print(jffpbruner.seqname+":"+str(jffpbruner.start)+':'+str(jffpbruner.end), "finished!")

    return probelist


def jfprobekmerfilter(jfpbkfruner):

    """
    :param jfpath: jellyfish bin path
    :param jfkmerfile: jellyfish kmer count file
    :param mer: int, kmer
    :param sequence: string, sequence for kmerscore count
    :param max: max kmer score
    :param min: min kmer score
    :return: list, kmerscore list
    """
    # jfpath = , jfkmerfile, mer, probe, maxk, mink

    probeinfo = jfpbkfruner.probe.split('\t')

    sequence = probeinfo[3]

    seqlen = len(sequence)

    jfpath = subprocesspath.subprocesspath(jfpbkfruner.jfpath)

    jfkmerfile = subprocesspath.subprocesspath(jfpbkfruner.jfkmerfile)

    jfquerycommand = ' '.join([jfpath, 'query', '-i', '-l', jfkmerfile])

    # print(jfquerycommand)

    kmerct = subprocess.Popen(jfquerycommand, shell=True, stdout=subprocess.PIPE,
                              stdin=subprocess.PIPE)


    mer = int(jfpbkfruner.mer)

    end = mer

    jfkmercount = list()

    keep = True

    while (end <= seqlen):

        start = end - mer

        subseq = sequence[start:end]+'\n'

        kmerct.stdin.write(subseq.encode('ascii'))

        kmerct.stdin.flush()

        lin = kmerct.stdout.readline().decode('utf-8').rstrip('\n')

        number = int(lin)

        # print(number)

        if number >= jfpbkfruner.maxk:

            keep = False

        if number <= jfpbkfruner.mink:

            keep = False

        jfkmercount.append(number)

        end += 1

    kmerct.stdin.close()

    kmerct.stdout.close()

    # kmerct.terminate()

    kmerct.wait()

    jfprobefileter = dict()

    jfprobefileter['chro'] = probeinfo[0]
    jfprobefileter['start'] = probeinfo[1]
    jfprobefileter['end'] = probeinfo[2]
    jfprobefileter['seq'] = probeinfo[3]
    jfprobefileter['keep'] = keep
    jfprobefileter['sumscore'] = sum(jfkmercount)

    return jfprobefileter


def jfseqNGSkmer(jfpath, jfkmerfile, mer, sequence, bfcount=False):

    """
    :param jfpath: jellyfish bin path
    :param jfkmerfile: jellyfish kmer count file
    :param mer: int, kmer
    :param sequence: string, sequence for kmerscore count
    :param bfcount:
    :return: list, kmerscore list
    """

    seqlen = len(sequence)

    jfpath = subprocesspath.subprocesspath(jfpath)

    jfkmerfile = subprocesspath.subprocesspath(jfkmerfile)

    jfquerycommand = ' '.join([jfpath, 'query', '-i', '-l', jfkmerfile])

    print(jfquerycommand)

    kmerct = subprocess.Popen(jfquerycommand, shell=True, stdout=subprocess.PIPE,
                              stdin=subprocess.PIPE)


    mer = int(mer)

    end = mer

    jfkmercount = list()

    while (end <= seqlen):

        start = end - mer

        subseq = sequence[start:end]+'\n'

        kmerct.stdin.write(subseq.encode('ascii'))

        kmerct.stdin.flush()

        lin = kmerct.stdout.readline().decode('utf-8').rstrip('\n')

        number = int(lin)

        jfkmercount.append(number)
        end += 1

    kmerct.stdin.close()

    kmerct.stdout.close()

    kmerct.wait()

    return jfkmercount


def jfngsscoer(jfngsscoer):

    """
    :param jfngsscoer: JFNGSScoer without score
    :return:  jfngsscoer with score
    """

    sequence = jfngsscoer.pyfasta[jfngsscoer.seqfullname][jfngsscoer.start:jfngsscoer.end]

    print("NGS scorer ", jfngsscoer.seqname+":"+str(jfngsscoer.start)+':'+str(jfngsscoer.end))

    jfngsscoer.score = jfseqNGSkmer(jfngsscoer.jfpath, jfngsscoer.jfkmerfile, jfngsscoer.mer, sequence)

    return jfngsscoer


def jfngsscoerlargegenome(jfngsscoer, tmppath):

    """
    :param jfngsscoer: JFNGSScoer without score
    :return:  jfngsscoer with score
    """

    sequence = jfngsscoer.pyfasta[jfngsscoer.seqfullname][jfngsscoer.start:jfngsscoer.end]

    print("NGS scorer ", jfngsscoer.seqname+":"+str(jfngsscoer.start)+':'+str(jfngsscoer.end))

    score = jfseqNGSkmer(jfngsscoer.jfpath, jfngsscoer.jfkmerfile, jfngsscoer.mer, sequence)

    tmpfile = jfngsscoer.seqname+'_'+str(jfngsscoer.start)+"_"+str(jfngsscoer.end)

    tmpfilename = os.path.join(tmppath, tmpfile)

    tmpio = open(tmpfilename, 'w')

    tmpio.write(' '.join(map(str,score)))
    #return jfngsscoer
    tmpio.close()

    print('save to ', tmpfilename)

class JFfpbruner():

    def __init__(self, jfpath, jfkmerfile, mer, pyfasta, seqname, start, end, pblength, maxkmerscore, step=4):

        self.jfpath = jfpath

        self.jfkmerfile = jfkmerfile

        self.mer = mer


        self.pblength = pblength

        self.maxkmerscore = maxkmerscore

        self.step = step

        self.pyfasta = pyfasta

        self.seqname = seqname

        self.start = start

        self.end =end


class JFpbkfruner():

    def __init__(self, jfpath, jfkmerfile, mer, probe, maxk, mink):

        self.jfpath = jfpath

        self.jfkmerfile = jfkmerfile

        self.mer = mer

        self.probe = probe

        self.maxk = maxk

        self.mink = mink


class JFNGSScoer():

    def __init__(self, jfpath, pyfasta, jfkmerfile, mer, seqfullname, start, end):
        self.jfpath = jfpath

        self.jfkmerfile = jfkmerfile

        self.mer = mer

        self.pyfasta = pyfasta

        self.seqfullname = seqfullname

        self.start = start

        self.end = end

        # convert name like '10 dna:chromosome chromosome:AGPv3:10:1:149632204:1' to '10'
        self.seqname = seqfullname.split()[0]

        self.score = list()


if __name__ == '__main__':
    import loadfa
#
    jfpath = '/Users/Forrest/Box Sync/Project/Chorus/bin/jellyfish/x86_64-Darwin/bin/jellyfish'

    print(jfversion(jfpath))

    jffile = '/Users/Forrest/Box Sync/Project/Chorus/Test/Testsampe/DM_404.fa_17mer.jf'

    sequence = loadfa.loadfa('/Users/Forrest/Box Sync/Project/Chorus/Test/DM_test.fa')

    jffbpruner = JFfpbruner(jfpath=jfpath, jfkmerfile=jffile, mer=17, sequence=sequence, pblength=45, maxkmerscore=33)

    pblist = kmerfilterprobe(jffbpruner)

    for pb in pblist:

        print(pb)
    #
    # jfscore = jfseqkmercount(jfpath=jfpath, jfkmerfile=jffile, mer=17, sequence='ATTGACAGGCGTTCAGGTAAGGAACTTGAGAA',
    #                          )
    #
    # print(jfscore, 0, 'to',10,sum(jfscore[0:10]))
#
#     outfile = 'testjf_33m.jf'
#
#     infile = '../testdata/chrM.fas'
#
#     jfcount(jfpath=jfpath, mer=33, threads=4, output=outfile, infile=infile, uppercount=1, size=10000000)
#
#     print(jfquery(jfpath=jfpath, countfile=outfile, sequence='ACAACATGGTTCCCGGTCGCGCATCCCCATGGA'))
#
#     print(jfquery(jfpath=jfpath, countfile=outfile, sequence='TGAGTAAAGCAGCAAATTTGTTTCACTAGCAAA'))

