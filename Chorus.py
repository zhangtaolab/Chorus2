import argparse
import sys
from Choruslib import bwa
from Choruslib import jellyfish
from Choruslib import prefilter, primer3_filter, spgenome
import os
from multiprocessing import Pool
from math import log
from pyfasta import Fasta


def main():

    args = check_options(get_options())

    genomesize = int(os.path.getsize(args.genome)/1e6)

    kmer = int(log(genomesize, 4)+1)

    if kmer < 17:

        kmer = 17

    #jellyfish par
    lowercount = 2

    #jellyfish par
    jfsize = '100M'

    # splite sequence longer than 10M
    spsize = 10000000

    step = args.step

    maxkmerscore = int(((args.length * args.homology / 100) - kmer) * args.ploidy/2 + 0.5 )

    jfpool = Pool(args.threads)

    # ?build kmerindex
    jfkmerfile = os.path.join(args.saved,(os.path.basename(args.genome)+'_'+str(kmer)+'mer.jf'))

    kmerbuild = True

    if os.path.isfile(jfkmerfile):

        if not args.docker:

            print("find:", jfkmerfile)

            kmmess = "Found kmerfile "+jfkmerfile+". Do you want rebuild it?  Press Y or N to continue:"

            print(kmmess)

            while True:

                char = getch()

                if char.lower() in ("y", "n"):

                    print(char)

                    if char == 'y':

                        kmerbuild = True

                    elif char == 'n':

                        kmerbuild = False

                    break


    # ?build bwa index
    bwaindexfile = os.path.basename(args.genome)

    bwatestindex = os.path.join(args.saved, bwaindexfile+'.sa')

    bwaindex = os.path.join(args.saved, bwaindexfile)

    bwabuild = True

    if os.path.isfile(bwatestindex):

        if not args.docker:

            print('find:', bwatestindex)

            bwamess = "Found bwa index file " + bwatestindex + ". Do you want rebuild it? Press Y or N to continue:"

            print(bwamess)

            while True:

                char = getch()

                if char.lower() in ("y", "n"):

                    print(char)

                    if char == 'y':

                        bwabuild = True

                    elif char == 'n':

                        bwabuild = False

                    break

    print("genomesize:",genomesize, "kmer:",kmer, "jfkmerfile:",
          jfkmerfile, "kmerbuild:", kmerbuild, "bwabuild:", bwabuild, "threads:", args.threads)

    # Build Jellyfish index
    if kmerbuild:

        jfcount = jellyfish.jfcount(jfpath=args.jellyfish, mer=kmer, infile=args.genome, output=jfkmerfile,
                                    threads=args.threads, lowercount=lowercount, size=jfsize)

        if jfcount:

            print("JellyFish Count finished ...")

        else:

            print("JellyFish Count Error!!!")

            sys.exit(1)

    else:

        print("Use ", jfkmerfile)
    # End build Jellyfish index

    if bwabuild:

        bwa.bwaindex(args.bwa, args.genome, args.saved)

        print("bwa index build finished ...")

    else:

        print("Use", bwatestindex)


    jffilteredprobe = list()

#####

    if genomesize < 1000:

        fastain = Fasta(args.input)

        jffpbrunerlist = list()

        for seqname in fastain.keys():

            chrlen = len(fastain[seqname])

            if chrlen < spsize:

                start = 0

                end = chrlen - 1

                jffpbruner = jellyfish.JFfpbruner(jfpath=args.jellyfish, jfkmerfile=jfkmerfile, mer=kmer,
                                                  pyfasta=fastain, seqname=seqname, pblength=args.length,
                                                  maxkmerscore=maxkmerscore, start=start,
                                                  end=end, step=step)

                jffpbrunerlist.append(jffpbruner)

            else:

                chrblock = int(chrlen/spsize) + 1

                for i in range(chrblock):

                    start = i * spsize

                    end = start + spsize - 1

                    if end >= chrlen:

                        end = chrlen - 1

                    jffpbruner = jellyfish.JFfpbruner(jfpath=args.jellyfish, jfkmerfile=jfkmerfile, mer=kmer,
                                                  pyfasta=fastain, seqname=seqname, pblength=args.length,
                                                  maxkmerscore=maxkmerscore, start=start,
                                                  end=end, step=step)

                    jffpbrunerlist.append(jffpbruner)

        jffinished = 0

        print(len(jffpbrunerlist))

        for curpblist in jfpool.imap_unordered(jellyfish.kmerfilterprobe, jffpbrunerlist):

            jffilteredprobe.extend(curpblist)

            jffinished += 1

            print("Jellyfish filter: ",jffinished,'/',len(jffpbrunerlist), sep='')

        jfpool.close()

        print('Jellyfish filter finished!!')

    else:

        ### split fa file when geome size greater than 1 Gb

        print("genome size > 1G")

        subFas = spgenome.spgenome(args.input, args.saved)



        for subFafile in subFas:
            print(subFafile)
            fastain = Fasta(subFafile)

            jffpbrunerlist = list()

            for seqname in fastain.keys():

                chrlen = len(fastain[seqname])

                if chrlen < spsize:

                    start = 0

                    end = chrlen - 1

                    jffpbruner = jellyfish.JFfpbruner(jfpath=args.jellyfish, jfkmerfile=jfkmerfile, mer=kmer,
                                                      pyfasta=fastain, seqname=seqname, pblength=args.length,
                                                      maxkmerscore=maxkmerscore, start=start,
                                                      end=end, step=step)

                    jffpbrunerlist.append(jffpbruner)

                else:

                    chrblock = int(chrlen / spsize) + 1

                    for i in range(chrblock):

                        start = i * spsize

                        end = start + spsize - 1

                        if end >= chrlen:
                            end = chrlen - 1

                        jffpbruner = jellyfish.JFfpbruner(jfpath=args.jellyfish, jfkmerfile=jfkmerfile, mer=kmer,
                                                          pyfasta=fastain, seqname=seqname, pblength=args.length,
                                                          maxkmerscore=maxkmerscore, start=start,
                                                          end=end, step=step)

                        jffpbrunerlist.append(jffpbruner)

            jffinished = 0

            print(len(jffpbrunerlist))

            for curpblist in jfpool.imap_unordered(jellyfish.kmerfilterprobe, jffpbrunerlist):
                jffilteredprobe.extend(curpblist)

                jffinished += 1

                print(subFafile + " Jellyfish filter: ", jffinished, '/', len(jffpbrunerlist), sep='')


        jfpool.close()

        print('Jellyfish filter finished!!')


    tmppbfa = os.path.join(args.saved, os.path.basename(args.input)+'_tmp_probe.fa')

    tmppbfaio = open(tmppbfa, 'w')

    seqnum = 0

    for tmppb in jffilteredprobe:

        print('>','seq',seqnum, sep='',file=tmppbfaio)

        print(tmppb,file=tmppbfaio)

        seqnum += 1

    tmppbfaio.close()

    del jffilteredprobe

    bwafiltedpb = bwa.bwafilter(bwabin=args.bwa, reffile=bwaindex, inputfile=tmppbfa, minas=args.length,
                                maxxs=int(args.length*args.homology/100), threadnumber=args.threads)

    # print(bwafiltedpb)

    tmpbwaftlist = os.path.join(args.saved, os.path.basename(args.input)+'.bed')

    alltmpbwaftlist = os.path.join(args.saved, os.path.basename(args.input)+'_all.bed')

    tmpbwaftlistio = open(tmpbwaftlist,'w')

    allbwaftlistio = open(alltmpbwaftlist,'w')

    seqlenfile = os.path.join(args.saved, os.path.basename(args.input)+'.len')

    seqlenio = open(seqlenfile,'w')

    seqlength = bwa.bwareflength(bwabin=args.bwa, reffile=bwaindex)

    for seqname in seqlength:

            print(seqname, seqlength[seqname], sep='\t', file=seqlenio)

    seqlenio.close()


    oligobefortmf = list()

    for pbtmp in bwafiltedpb:

        # print(pbtmp, file=tmpbwaftlistio)
        nowpbcounter = dict()

        nowpbcounter['seq'] = pbtmp

        nowpbcounter['dTm'] = args.dtm

        nowpbcounter['rprimer'] = args.primer

        oligobefortmf.append(nowpbcounter)

    keepedprobe = list()

    ctedpb = 0

    oligobefortmflen = len(oligobefortmf)

    print("oligobefortmflen:",oligobefortmflen)

    pbftpool = Pool()

    for (pb, keep) in pbftpool.imap_unordered(probefilter, oligobefortmf):

        if keep:

            keepedprobe.append(pb)
                # print(pb, file=tmpbwaftlistio)
        ctedpb += 1

        if ctedpb % 10000 == 0:

            print(ctedpb,'/',oligobefortmflen)

    pbdictbychr = dict()

    pbftpool.close()

    for pb in keepedprobe:

        seq, chro, start = pb.split('\t')

        start = int(start)

        if chro in pbdictbychr:

            pbdictbychr[chro][start] = seq

        else:

            pbdictbychr[chro] = dict()

            pbdictbychr[chro][start] = seq

    lenrprimer = len(args.primer)

    if lenrprimer == 0:

            lenrprimer = 5

    slidwindow = lenrprimer+args.length

    for chro in pbdictbychr:

        startn = 0

        for startnow in sorted(pbdictbychr[chro]):

            endnow = startnow + args.length - 1

            print(chro, startnow, endnow, pbdictbychr[chro][startnow],file=allbwaftlistio,sep='\t')

            if startnow > startn+slidwindow:
                    #startn = startnow+slidwindow
                startn = startnow

                print(chro, startnow, endnow, pbdictbychr[chro][startnow], file=tmpbwaftlistio, sep='\t')


    tmpbwaftlistio.close()

    allbwaftlistio.close()

    print("Job finshed!!")

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

def check_options(parser):

    args = parser.parse_args()

    # Start check bwa
    if args.bwa:

        if not os.path.exists(args.bwa):

            print("Can not locate bwa, please input full path of bwa\n")

            parser.print_help()

            sys.exit(1)

        bwaversion = bwa.bwaversion(args.bwa)

        if bwaversion == 'None':

            print("Can not locate bwa, please input full path of bwa\n")

            parser.print_help()

            sys.exit(1)

    else:

        bwapath = which('bwa')

        if bwapath:

            bwaversion = bwa.bwaversion(bwapath[0])

            if bwaversion == 'None':

                print("Can not locate bwa, please input full path of bwa\n")

                parser.print_help()

                sys.exit(1)

            else:

                args.bwa = bwapath[0]

        else:

            print("Can not locate bwa, please input full path of bwa\n")

            parser.print_help()

            sys.exit(1)

    # End check bwa

    # Start check jellyfish
    if args.jellyfish:

        if not os.path.exists(args.jellyfish):

            print("Can not locate jellyfish, please input full path of jellyfish\n")

            parser.print_help()

            sys.exit(1)

        jellyfishversion = jellyfish.jfversion(args.jellyfish)

        if jellyfishversion == 'None':

            print("Can not locate jellyfish, please input full path of jellyfish\n")

            parser.print_help()

            sys.exit(1)

    else:

        jellyfishpath = which('jellyfish')

        if jellyfishpath:

            jellyfishversion = jellyfish.jfversion(jellyfishpath[0])

            if jellyfishversion == 'None':

                print("Can not locate jellyfish, please input full path of jellyfish\n")

                parser.print_help()

                sys.exit(1)

            else:

                args.jellyfish = jellyfishpath[0]

        else:

            print("Can not locate jellyfish, please input full path of jellyfish\n")

            parser.print_help()

            sys.exit(1)
    # End check jellyfish

    if not os.path.exists(args.genome):

        print("Can not locate genome file, please input genome file.\n")

        parser.print_help()

        sys.exit(1)

    if not os.path.exists(args.input):

        print("Can not locate input file, please input input file.\n")

        parser.print_help()

        sys.exit(1)

    if args.homology > 100 or args.homology < 50:

        print("Please set homology between 50 and 100.\n")

        parser.print_help()

        sys.exit(1)

    if args.dtm > 37 or args.dtm < 0:

        print("Please set dtm between 0 and 37.\n")

        parser.print_help()

        sys.exit(1)

    if args.step < 1:

        print("step must >1.\n")

        parser.print_help()

        sys.exit(1)

    # Start check saved folder
    if os.path.exists(args.saved):

        if not args.docker:

            print(args.saved, "exists. Everything in this folder will be remove. Press Y or N to continue: ")

            while True:

                char = getch()

                if char.lower() in ("y", "n"):

                    print(char)

                    if char == 'n':

                        sys.exit(1)

                    break

    else:

        os.mkdir(args.saved)
    # End check saved folder

    # Print Checked information
    print("#"*40)

    print("bwa version:", args.bwa, bwaversion)

    print("jellyfish version:", args.jellyfish, jellyfishversion)

    print("genome file:", args.genome)

    print("input file:", args.input)

    print("5\' labeled R primer:", args.primer)

    print("result output folder:",  os.path.realpath(args.saved))

    print("threads number:", args.threads)

    print("homology:", args.homology)

    print("dtm:", args.dtm)

    print("#"*40)

    return args


def getch():
    """
    For yes/no choice
    """
    import sys, tty, termios
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        ch = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return ch


def which(filename):
    """docstring for which"""
    locations = os.environ.get("PATH").split(os.pathsep)
    candidates = []
    for location in locations:
        candidate = os.path.join(location, filename)
        if os.path.isfile(candidate):
            candidates.append(candidate)
    return candidates


def get_options():

    parser = argparse.ArgumentParser(description="Chorus Software for Oligo FISH probe design", prog='Chorus',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="Example:\n"
                                            "  python3 Chorus.py -i TAIR10_chr_all.fas -g TAIR10_chr_all.fas -t 4 \\ \n"
                                            "                    -j /opt/software/jellyfish/bin/jellyfish -b /opt/software/bwa/bwa -s sample"
                                     )

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('-j', '--jellyfish', dest='jellyfish', help='The path where Jellyfish software installed')

    parser.add_argument('-b', '--bwa', dest='bwa', help='The path where BWA software installed')

    parser.add_argument('-g', '--genome', dest='genome', help='Fasta format genome file, should include all sequences from genome', required=True)

    parser.add_argument('-i', '--input', dest='input', help='Fasta format input file, can be whole genome, a chromosome or one region from genome', required=True)

    parser.add_argument('-s', '--save', dest='saved', help='The output folder for saving results', default='probes')

    parser.add_argument('-p', '--primer', dest='primer', help='A specific 5\' labeled R primer for PCR reaction. For example: CGTGGTCGCGTCTCA. (Default is none)', default='')

    parser.add_argument('-t', '--threads', dest='threads', help='Number of threads or CPUs to use. (Default: 1)',
                        default=1, type=int)

    parser.add_argument('-l', '--length', dest='length', help='The probe length. (Default: 45)', default=45, type=int)

    parser.add_argument('--homology', dest='homology', help='The maximum homology(%%) between target sequence and probe, range from 50 to 100. (Default: 75)', default=75, type=float)

    parser.add_argument('-d', '--dtm', dest='dtm', help='The minimum value of dTm (hybrid Tm - hairpin Tm), range from 0 to 37. (Default: 10)', default=10, type=float)

    parser.add_argument('--step', dest='step', help='The step length for k-mer searching in a sliding window, step length>=1. (Default: 5)', default=5, type=int)

    parser.add_argument('--docker', help='Only used in Docker version of Chorus', default=False)

    parser.add_argument('--ploidy', dest='ploidy', default=2, type=int, help='The ploidy of the given genome (test version). (Default: 2)')
    # parser.parse_args(['--version'])
    # args = parser.parse_args()

    return parser

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)
