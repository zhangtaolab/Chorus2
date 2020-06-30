import sys
import argparse
from Choruslib import bwa, bcftools, probecompare, bamdepth, jellyfish, subprocesspath, revcom
import os
from pybedtools import BedTool
from multiprocessing import Pool
import pandas as pd

def main():

    args = check_options(get_options())

    # jellyfish par
    jfsize = '100M'

    # ?build bwa index
    bwaindexfile = os.path.basename(args.genome)

    tmpfolder = args.tmp

    bwatestindex = os.path.join(tmpfolder, bwaindexfile+'.sa')

    bwaindex = os.path.join(tmpfolder, bwaindexfile)



    """
        bwabuild = True
    
        if os.path.isfile(bwatestindex):
    
            bwabuild = False
    
        if bwabuild:
    
            # build bwa index
            bwa.bwaindex(args.bwa, args.genome, tmpfolder)
    
            print("bwa index build finished ...")
    
        else:
    
            print("Use", bwatestindex)
    """

    sampleinfor = dict()

    names = args.names.split(',')

    reads1 = args.reads1.split(',')

    reads2 = args.reads2.split(',')

    cnsfile = os.path.join(args.saved,'_'.join(names)+'_cns_probe.csv')

    print(cnsfile)

    cnsio = open(cnsfile, 'w')

    for i in range(len(names)):

        name = names[i]

        read1 = reads1[i]

        read2 = reads2[i]

        bamfile = os.path.join(tmpfolder,name+'.bam')

        bcffile = os.path.join(tmpfolder,name+'.bcf')

        jffile = os.path.join(tmpfolder,name+'.jf')

        cnsprobe = os.path.join(args.saved, name + '_probe.txt')

        # new add indel
        indelNprobe = os.path.join(args.saved, name + '_indel_probe.txt')

        mindepth = os.path.join(tmpfolder, name + '_mindepth.bed')

        if name in sampleinfor:

            print("error same name:", name)

        else:

            sampleinfor[name] = dict()

            sampleinfor[name]['read1'] = read1

            sampleinfor[name]['read2'] = read2

            sampleinfor[name]['bamfile'] = bamfile

            sampleinfor[name]['bcffile'] = bcffile

            sampleinfor[name]['jffile'] = jffile

            # sampleinfor[name]['kmerscore'] = kmerscore
            #
            # sampleinfor[name]['kmerscoreio'] = open(kmerscore, 'w')

            sampleinfor[name]['cnsprobe'] = cnsprobe

            sampleinfor[name]['cnsprobeio'] = open(cnsprobe, 'w')

            # new add indel
            sampleinfor[name]['indelNprobelist'] = list()
            sampleinfor[name]['indelNprobeio'] = open(indelNprobe, 'w')

            sampleinfor[name]['mindepth'] = mindepth
            """
            # run bwa mem
            bwa.bwamem_paired(bwabin=args.bwa,
                          samtoolsbin =args.samtools,
                          reffile=bwaindex,
                          outfile=bamfile,
                          inputfile1=read1,
                          inputfile2=read2,
                          samplename=name,
                          threadnumber = args.threads
                         )
            
            print("bwa mem", name, 'finished')
            
            
            # get min depth bed file
            bamdepth.bamdepthtobed(bamfile=bamfile, outbed=mindepth, mindepth=args.mindepth, minlength=200)

            print(mindepth, 'done')

            # generate bcf file from bam file
            bcftools.bamtobcf(bcfbin=args.bcftools,
                     reffile=bwaindex,
                     bamfile=bamfile,
                     outbcf=bcffile)

            print(bcffile, "done")


            # generate jf file

            jellyfish.makegenerator(filenames=[read1, read2], type='gz', generators='generators')

            jellyfish.jfgeneratorscount(jfpath=args.jellyfish, mer=args.length, output=jffile,
                                        generators='generators', threads=args.threads, size='100M')

            print(jffile, "done")
            """

    probe = BedTool(args.probe).sort()

    for name in sampleinfor:

        nowprobe = BedTool(sampleinfor[name]['mindepth']).sort()

        probe = probe.intersect(nowprobe, wa=True, u=True)

    # cnsprobe



    for name in sampleinfor:

        bcfpool = Pool(args.threads)

        bcfrunerlist = list()

        consensusprobelist = list()

        for i in probe:

            probestr = str(i).rstrip()

            bcfconsensusruner = bcftools.BcfConsensusRuner(probestr=probestr, bcftoolspath=args.bcftools,
                                                           bcffile=sampleinfor[name]['bcffile'],
                                                           sample=name
                                                           )

            bcfrunerlist.append(bcfconsensusruner)
            # consensusprobe = bcftools.probestrtoconsensus(bcfconsensusruner)
            #
            # print(probestr, consensusprobe, sep='\t')

        reslist = list()

        for res in bcfpool.imap_unordered(bcftools.probestrtoconsensus, bcfrunerlist):

            # print(res['probestr'], name, res['consensusprobe'], sep='\t', file=sampleinfor[name]['cnsprobeio'])

            if len(res['consensusprobe']) != args.length:

                sampleinfor[name]['indelNprobelist'].append(res)

            elif 'N' in res['consensusprobe']:

                continue

            else:

                consensusprobelist.append(res['consensusprobe'])
                # consensusprobelist.append(res)
                reslist.append(res)


        bcfpool.close()

        consensusprobekmerscore = jellyfish.jfquerylist(jfkmerfile=sampleinfor[name]['jffile'],
                                                        jfpath=args.jellyfish,
                                                        seqlist=consensusprobelist)

        kmerscoredict = dict()

        kmerscorelist = list()

        for score in consensusprobekmerscore:

            # print(score, file=sampleinfor[name]['kmerscoreio'])
            (subseq, kmerscore) = score.split(',')

            if 'N' not in subseq:

                kmerscoredict[subseq] = int(kmerscore)

                kmerscorelist.append(int(kmerscore))

        maxkmer = pd.Series(kmerscorelist).quantile(0.9)

        minkmer = args.minkmer

        for consensusprobe in reslist:

            probestr = consensusprobe['probestr']

            consensusprobeseq = consensusprobe['consensusprobe']

            if consensusprobeseq in kmerscoredict:

                if kmerscoredict[consensusprobeseq] <= maxkmer:

                    if kmerscoredict[consensusprobeseq] >= minkmer:

                        print(probestr, consensusprobeseq, kmerscoredict[consensusprobeseq], sep='\t',
                              file=sampleinfor[name]['cnsprobeio'])

    for name in sampleinfor:

        sampleinfor[name]['cnsprobeio'].close()
        # sampleinfor[name]['kmerscoreio'].close()
    # print(sampleinfor)

        for res in sampleinfor[name]['indelNprobelist']:

            print(res['probestr'], name, res['consensusprobe'], sep='\t', file=sampleinfor[name]['indelNprobeio'])

        sampleinfor[name]['indelNprobeio'].close()

    probdict = dict()

    for name in sampleinfor:

        with open(sampleinfor[name]['cnsprobe']) as inio:

            for infor in inio:

                infor = infor.rstrip()

                inforlist = infor.split('\t')

                orgprb = inforlist[3]

                if orgprb in probdict:

                    probdict[orgprb][name] = infor

                else:

                    probdict[orgprb] = dict()

                    probdict[orgprb][name] = infor

    print('chrom', 'start', 'end', 'refseq', ','.join(sampleinfor), 'consensusprobe', 'consensusscore', 'consensussite',
          'consensusdiff', sep=',', file=cnsio)

    for orgprb in probdict:

        sharecount = len(probdict[orgprb])

        values_view = probdict[orgprb].values()
        value_iterator = iter(values_view)
        first_value = next(value_iterator).split('\t')

        outinfo = first_value[0:3]

        if len(sampleinfor) == sharecount:
            #         print(sampleinfor, sharecount)
            # print(orgprb, len(probdict[orgprb]))
            probelist = list()
            namelist = list()
            namelist.append('refseq')
            probelist.append(orgprb)

            for name in sampleinfor:

                infor = probdict[orgprb][name].split('\t')

                speciesprobe = infor[-2]

                namelist.append(name)
                if len(speciesprobe) == len(orgprb):
                    probelist.append(speciesprobe)

            if len(namelist) == len(probelist):
                #             print(namelist, probelist)

                res = probecompare.getconsensusprobe(probelist)
                outinfo.extend(probelist)
                print(','.join(outinfo), res['consensusprobe'], res['consensusscore'], res['consensussite'],
                      res['consensusdiff'], sep=',', file=cnsio)

    cnsio.close()

    print("finished")

def check_options(parser):

    args = parser.parse_args()

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

        bwapath = subprocesspath.which('bwa')

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

        jellyfishpath = subprocesspath.which('jellyfish')

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

    if args.samtools:

        if not os.path.exists(args.samtools):

            print("Can not locate samtools, please input full path of samtools\n")

            parser.print_help()

            sys.exit(1)

        samtoolsversion = bwa.samtoolsversion(args.samtools)

        if samtoolsversion == 'None':

            print("Can not locate samtools, please input full path of samtools\n")

            parser.print_help()

            sys.exit(1)

    else:

        samtoolspath = subprocesspath.which('samtools')

        if samtoolspath:

            samtoolsversion = bwa.samtoolsversion(samtoolspath[0])

            if samtoolsversion == 'None':

                print("Can not locate samtools, please input full path of samtools\n")

                parser.print_help()

                sys.exit(1)

            else:

                args.samtools = samtoolspath[0]

        else:

            print("Can not locate samtools, please input full path of samtools\n")

            parser.print_help()

            sys.exit(1)

    if args.bcftools:

        if not os.path.exists(args.bcftools):

            print("Can not locate bcftools, please input full path of bcftools\n")

            parser.print_help()

            sys.exit(1)

        bcftoolsversion = bcftools.bcftoolsversion(args.bcftools)

        if bcftoolsversion == 'None':

            print("Can not locate bcftools, please input full path of bcftools\n")

            parser.print_help()

            sys.exit(1)

    else:

        bcftoolspath = subprocesspath.which('bcftools')

        if bcftoolspath:

            bcftoolsversion = bcftools.bcftoolsversion(bcftoolspath[0])

            if bcftoolsversion == 'None':

                print("Can not locate bcftools, please input full path of bcftools\n")

                parser.print_help()

                sys.exit(1)

            else:

                args.bcftools = bcftoolspath[0]

        else:

            print("Can not locate bcftools, please input full path of bcftools\n")

            parser.print_help()

            sys.exit(1)

    if not os.path.exists(args.genome):

        print("Can not locate genome file, please input genome file.\n")

        parser.print_help()

        sys.exit(1)

    if not os.path.exists(args.probe):

        print("Can not locate probe file, please input probe file which generate by Chorus2.\n")

        parser.print_help()

        sys.exit(1)

    names = args.names.split(',')

    reads1 = args.reads1.split(',')

    reads2 = args.reads2.split(',')

    if len(names) != len(reads1):

        print("There are %s species (%s) but %s reads1 files (%s)" % (len(names), names, len(reads1), reads1))

    if len(names) != len(reads2):
        print("There are %s species (%s) but %s reads2 files (%s)" % (len(names), names, len(reads2), reads2))


    for read1 in reads1:
        # print(read1)
        if not os.path.exists(read1):

            print("Can not locate %s file.\n", read1)

            parser.print_help()

            sys.exit(1)

    for read2 in reads2:

        if not os.path.exists(read2):
            print("Can not locate %s file.\n", read2)

            parser.print_help()

            sys.exit(1)

    if os.path.exists(args.saved):
            args.saved = os.path.realpath(args.saved)

            print("The output folder ",args.saved, " already exists.")
            print('''
                Press Y to use it for output files. Everything in this folder will be removed.
                Press N and set -s/--saved to a different folder:
                ''')

            while True:

                char = getch()

                if char.lower() in ("y", "n"):

                    print(char)

                    if char.lower() == 'n':

                        sys.exit(1)

                    break

    else:

        args.saved = os.path.realpath(args.saved)

        os.mkdir(args.saved)

    if args.tmp:

        tmpfolder = os.path.abspath(args.tmp)

    else:

        tmpfolder = os.path.abspath(os.path.join(args.saved, './tmp'))

    print("tmp",args.tmp, tmpfolder)

    args.tmp = tmpfolder

    if os.path.exists(tmpfolder):

            tmpfolder = os.path.realpath(tmpfolder)

            print("The output folder ",tmpfolder, " already exists.")
            print('''
                Press Y to use it for output files. Everything in this folder will be removed.
                Press N and set --tmp to a different folder:
                ''')

            while True:

                char = getch()

                if char.lower() in ("y", "n"):

                    print(char)

                    if char.lower() == 'n':

                        sys.exit(1)

                    break

    else:

        args.tmp = tmpfolder

        os.mkdir(tmpfolder)


    print("#"*40)

    print("bwa version:", args.bwa, bwaversion)

    print("jellyfish version:", args.jellyfish, jellyfishversion)

    print("samtools version:", args.samtools, samtoolsversion)

    print("bcftools version:", args.bcftools, bcftoolsversion)

    print("genome file:", args.genome)

    print("probe file:", args.probe)

    print("out put folder:", args.saved)

    for i in range(len(names)):

        print(names[i], reads1[i], reads2[i])

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



def get_options():
    parser = argparse.ArgumentParser(description="Oligo FISH probe design for no reference genome.", prog='ChorusNoRef',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="Example:\n"

                                     )

    parser.add_argument("--version", action='version', version='%(prog)s 2.1' )

    parser.add_argument('-j', '--jellyfish', dest='jellyfish', help='The path where Jellyfish software installed')

    parser.add_argument('-b', '--bwa', dest='bwa', help='The path where BWA software installed')

    parser.add_argument('-c', '--bcftools', dest='bcftools', help='The path where bcftools software installed')

    parser.add_argument('-m', '--samtools', dest='samtools', help='The path where samtools software installed')

    parser.add_argument('-g', '--genome', dest='genome',
                        help='Fasta format genome file, should include all sequences from genome', required=True)

    parser.add_argument('-s', '--save', dest='saved', help='The output folder for saving results', default='noRefprobes')

    parser.add_argument('--tmp', dest='tmp', help='The temporary fold for processing')

    parser.add_argument('-p', '--probe', dest='probe',
                        help='Original probe design by Chorus2',
                        required=True)

    parser.add_argument('-r1','--reads1',dest='reads1', required=True,
                        help='read1 of species, example: For one Species only: species_R1.fq; for more than one species: sepecies1_R1.fq,sepecies2_R1.fq ')

    parser.add_argument('-r2', '--reads2', dest='reads2', required=True,
                        help='read1 of species, example: For one Species only: species_R2.fq; for more than one species: sepecies1_R2.fq,sepecies2_R2.fq ')

    parser.add_argument('-n', '--names', dest='names', required=True,
                        help='species name(s), the order must same as r1, r2 ')

    parser.add_argument('-t', '--threads', dest='threads', help='Number of threads or CPUs to use. (Default: 1)',
                        default=1, type=int)

    parser.add_argument('--minkmer', default=3, type=int, dest='minkmer', help="probe min count for illumina reads")

    parser.add_argument('-l', '--length', dest='length', help='The probe length. (Default: 45)', default=45, type=int)

    parser.add_argument('-d', '--mindepth', dest='mindepth', help='Mini depth covered by illumina sequence. (Default: 3)', default=3, type=int)

    return parser

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)


# python /mnt/e/Github/Chorus2/ChorusNoRef.py -g /mnt/e/Data/Ref/Potato/potato_dm_v404_all_pm_un.fasta -p /mnt/e/Data/Solanum/Oligo/potato_dm_v404/potato_dm_v404_all_pm_un_jffilted_pq2575.bed -r1 /mnt/e/Data/Solanum/Solanum_jamesii/SRR5349574_1.fastq.gz,/mnt/e/Data/Solanum/Solanum_etuberosum/SRR5349573_1.fastq.gz -r2 /mnt/e/Data/Solanum/Solanum_jamesii/SRR5349574_2.fastq.gz,/mnt/e/Data/Solanum/Solanum_etuberosum/SRR5349573_2.fastq.gz -n jamesii,etuberosum