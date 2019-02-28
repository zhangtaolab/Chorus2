import argparse
import sys
from Choruslib import jellyfish
import os
from multiprocessing import Pool, Process
from pyfasta import Fasta
import pyBigWig
import math

def main():

    args = check_options(get_options())

    # jfgeneratorscount(jfpath, mer, output, generators,threads=1,  size='100M'):

    # make generators

    print(args.input)

    jellyfish.makegenerator(filenames=args.input.split(','), type=args.gzip, generators='generators')

    jfkmerfile = args.output+'.jf'

    bwfile = args.output+'.bw'

    outfilename = args.output

    jellyfish.jfgeneratorscount(jfpath=args.jellyfish, mer=args.kmer, output=jfkmerfile,
                                generators='generators',threads=args.threads,  size='100M')

    spsize = 10000000

    fastain = Fasta(args.genome)

    bw = pyBigWig.open(bwfile, "w")

    seqlenth = dict()
    seqname = dict()

    genomesize = 0

    for chrom in sorted(fastain.keys()):
        infor = chrom.split()
        seqlenth[infor[0]] = len(fastain[chrom])
        seqname[infor[0]] = chrom
        genomesize += seqlenth[infor[0]]

    print("Genome Size: %s" % genomesize)
    bw.addHeader(list(seqlenth.items()))

    jfscoerlist = list()


    for seqfullname in sorted(fastain.keys()):

        infor = seqfullname.split()

        chrlen = len(fastain[seqfullname])

        if chrlen < spsize:

            start = 0

            end = chrlen - 1

            jfscoer = jellyfish.JFNGSScoer(jfpath=args.jellyfish, jfkmerfile=jfkmerfile, mer=args.kmer,
                                           start=start, end=end, seqfullname=seqfullname, pyfasta=fastain)

            jfscoerlist.append(jfscoer)

        else:

            chrblock = int(chrlen / spsize) + 1

            for i in range(chrblock):

                start = i * spsize

                end = start + spsize - 1

                if i > 0:
                    start = start - args.kmer + 1

                if end >= chrlen:
                    end = chrlen - 1

                jfscoer = jellyfish.JFNGSScoer(jfpath=args.jellyfish, jfkmerfile=jfkmerfile, mer=args.kmer,
                                               start=start, end=end, seqfullname=seqfullname, pyfasta=fastain)
                jfscoerlist.append(jfscoer)


    tmppath = os.path.dirname(args.output)

    jfsllength = int(len(jfscoerlist) / args.threads + 1)

    for jt in range(jfsllength):

        if jt == jfsllength:

            nowlist = jfscoerlist[jt * args.threads:]

        else:

            nowlist = jfscoerlist[(jt * args.threads):((jt + 1) * args.threads)]

        processes = list()

        for jfscoer in nowlist:

                p = Process(target=jellyfish.jfngsscoerlargegenome, args=(jfscoer,tmppath))

                processes.append(p)

        for p in processes:

                p.start()

        for p in processes:

                p.join()
            # jfngsscoerlargegenome

        for jfscoer in nowlist:

                tmpfile = jfscoer.seqname + '_' + str(jfscoer.start) + "_" + str(jfscoer.end)

                tmpfilename = os.path.join(tmppath, tmpfile)

                score = list()

                with open(tmpfilename) as inio:

                    for i in inio:

                        score = i.rstrip().split()

                bw.addEntries(jfscoer.seqname, jfscoer.start, values=list(map(float, score)), span=1, step=1)

                print(jfscoer.seqname, jfscoer.start, 'OK')

                inio.close()

                os.remove(tmpfilename)

    bw.close()


    bwforcount = pyBigWig.open(bwfile)



    outio = open(outfilename, 'w')

    with open(args.probe) as inio:

        for i in inio:

            #(chrom, start, end, seq) = i.rstrip().split()
            probeloc = i.rstrip().split()

            chrom = probeloc[0]

            start = probeloc[1]

            end = probeloc[2]

            seq = probeloc[3]

            score = sum(bwforcount.values(chrom, int(start) - 1, int(end) - args.kmer))

            if math.isnan(score):

                score = 0

            print(chrom, start, end, seq, int(score), '+', sep='\t', file=outio)

    outio.close()

    print("finished!")

def check_options(parser):

    args = parser.parse_args()

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

    if not os.path.exists(args.probe):

        print("Can not locate probe file, please input genome file.\n")

        parser.print_help()

        sys.exit(1)


    if args.input:

        inputfiles = args.input.split(',')

        for inputfile in inputfiles:

            if not os.path.exists(inputfile):

                print("Can not locate %s file.\n" % inputfile)

                parser.print_help()

                sys.exit(1)

    else:

        print("Can not locate input file, please input input file.\n")

        parser.print_help()

        sys.exit(1)

    return args

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

    parser = argparse.ArgumentParser(description="Chorus Software for counting Oligo FISH probe k-mer score using NGS data", prog='ChorusNGSfilter',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="Example:\n"
                                            "  python3 ChorusNGSfilter.py -i 1.fq.gz,2.fq.gz -z gz -t 4 -g TAIR10_chr_all.fas \\ \n"
                                            "                             -j /opt/software/jellyfish/bin/jellyfish -p probe.bed -o output.bed"
                                     )

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('-j', '--jellyfish', dest='jellyfish', help='The path where Jellyfish software installed')

    parser.add_argument('-g', '--genome', dest='genome', help='Fasta format genome file, should include all sequences from genome', required=True)

    parser.add_argument('-i', '--input', dest='input',
                        help='Fastq format input files contain reads from whole genome shotgun sequencing, files can be gzipped.'
                             ' Multiple files separate with \",\". For example: 1.fq.gz,2.fq.gz ', required=True, type=str)

    parser.add_argument('-z', '--gzipped', dest='gzip', help='Input fastq file is gzipped(gz) or uncompressed(text). (Default: gz)', choices=('gz', 'text'), default='gz', required=True)

    # parser.add_argument('-s', '--save', dest='saved', help='result saved folder', default='probes')

    parser.add_argument('-t', '--threads', dest='threads', help='Number of threads or CPUs to use. (Default: 1)',
                        default=1, type=int)

    parser.add_argument('-k', '--kmer', dest='kmer', help='Length of k-mer used for counting k-mers in input fastq files. (Default: 17)', default=17, type=int)

    parser.add_argument('-p', '--probe', dest='probe', help='The bed format probe file generated by Chorus')

    parser.add_argument('-o', '--output', dest='output', help='Output bed format probe file with k-mer score. (Default: output.bed)', default='output.bed')

    # args = parser.parse_args()

    return parser

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)
