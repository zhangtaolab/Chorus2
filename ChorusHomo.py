import os
import sys
import argparse
from subprocess import call
import re

from Choruslib import jellyfish
from Choruslib import bwa
from Choruslib import subprocesspath


def main():
    # get args
    args = check_options(get_options())
    print(args)

    # init input and output files
    related_genome = args.genome
    target = args.target

    probe_file = args.input
    probe_fasta = bed_to_fa(probe_file, args.saved)
    homo_probe = make_saved_file(args.saved, probe_file, "homo.csv")

    sys.stderr.write("input probe {}\n".format(probe_file))

    # Check related genome bwa index or build it
    if not check_bwa_index(related_genome, args.bwa, args.saved):
        sys.stderr.write("Failed to get index for ", related_genome, "\n")
        sys.exit(0)

    indexed_genome = os.path.join(args.saved, os.path.basename(related_genome))

    # generate fai file for each genome
    make_fai(indexed_genome)
    make_fai(target)

    tmp_sam_file = os.path.join(args.saved, 'tmp_align.sam')

    # align probes to genome
    bwafiltedpb = bwa_mem(
        bwabin=args.bwa,
        reffile=indexed_genome,
        inputfile=probe_fasta,
        outfile=tmp_sam_file,
        threadnumber=args.threads)
    #
    hom_probe_dict = dict()
    for op in bwafiltedpb:
        idx, probeseq, query_chr, query_st, query_ed, fake_num,  seqname, start, end, identity = op.strip().split(',')
        head = query_chr + query_st
        hom_probe_dict[head] = op

# add homology information after each oligos
    with open(homo_probe, "w") as o:
        with open(probe_file, "r") as i:
            for line in i:
                arr = line.strip().split()
                hh = arr[0] + arr[1]
                if hh in hom_probe_dict:
                    save_probe([hom_probe_dict[hh]], o)

    sys.stderr.write("record homologous oligos\n")


def make_fai(genome):
    cmd = 'samtools faidx ' + genome
    call(cmd, shell=True)


def bwa_mem(bwabin, reffile, inputfile, outfile, threadnumber=1):

    pat = re.compile('^@')

    # bwabin = subprocesspath.subprocesspath(bwabin)

    # reffile = subprocesspath.subprocesspath(reffile)

    # inputfile = subprocesspath.subprocesspath(inputfile)

    # bwacmd = ' '.join([bwabin, 'mem', '-O',' 0',' -B',' 0',' -E',' 0',' -k',' 5', '-t',str(threadnumber), reffile, inputfile])

    # print(bwacmd)

    bwa.bwaalign(bwabin, reffile, inputfile, outfile, threadnumber)

#    aspat = re.compile('AS:i:(\d.)')
#
#    xspat = re.compile('XS:i:(\d.)')

    res = list()
    idx=0

    tmp_sam_in = open(outfile, 'r')

    for lin in tmp_sam_in.readlines():
        # print("before decode",lin)
        lin = lin.rstrip('\n')
        # print("after decode", lin)
        if not re.search(pat, lin):
            
            infor = lin.split('\t')
            idx = idx + 1
            query_name = infor[0]
            query_chr,query_st,query_ed=query_name.split('_')
            seqname = infor[2]

            start = infor[3]

            probeseq = infor[9]

            md = re.split(':',infor[12])[-1]

            aln_matches = sum([int(item) for item in re.split('[ACTG^]', md) if not item == ''])
            aln_mismatches = sum([len(item) for item in re.split('[\d+^]', md) if not item == ''])

            if (aln_matches + aln_mismatches) == 0:

                identity = 0

            else:

                identity = aln_matches / (aln_matches + aln_mismatches)

#            asmatch = re.search(aspat, lin)
#
#            xsmatch = re.search(xspat, lin)
#
#            if asmatch:
#
#                asscore = int(asmatch.group(1))
#
#            else:
#
#                continue
#
#            if xsmatch:
#
#                xsscore = int(xsmatch.group(1))
#
#            else:
#
#                continue
#
#            if (asscore >= minas) & (xsscore < maxxs):

            end = str(int(start) + aln_matches + aln_mismatches -1)
            res.append(','.join([str(idx), probeseq, query_chr, query_st, query_ed, '0.99',
                                 seqname, start, end, str(f'{identity:.2f}')]))

    tmp_sam_in.close()

    os.remove(outfile)

    return res


def make_saved_file(saved_path, file, extention):
    ''' add full path to file '''
    name = os.path.basename(file)
    base = os.path.splitext(name)[0]
    saved_file = "".join((saved_path, "/", base + "." + extention))
    return saved_file


def check_bwa_index(related_genome, bwabin, saved_path):
    ''' check if bwa index exist '''
    index = related_genome + ".sa"
    if os.path.isfile(index):
        sys.stderr.write("index already existed\n")
        return True
    else:
        sys.stderr.write("make index file for genome\n")
        bwa.bwaindex(bwabin, related_genome, saved_path)
        return True

    return False


def bed_to_fa(probe_file, saved_path):
    ''' convert bed to fasta, header is the position of the seq.
    All other information in the bed file is lost'''
    fasta = make_saved_file(saved_path, probe_file, "fa")
    with open(fasta, "w") as o:
        with open(probe_file, "r") as i:
            for line in i:
                infor = line.strip().split()
                chrom = infor[0]
                st = infor[1]
                ed = infor[2]
                seq = infor[3]
                header = "_".join((chrom, st, ed))
                o.write(">" + header + "\n" + seq + "\n")

    return fasta


def save_probe(one_oligo, output_handle):
    ''' save oligo to file '''
    output_handle.write("\t".join(one_oligo))
    output_handle.write("\n")


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

            print(
                "Can not locate jellyfish, please input full path of jellyfish\n"
            )

            parser.print_help()

            sys.exit(1)

        jellyfishversion = jellyfish.jfversion(args.jellyfish)

        if jellyfishversion == 'None':

            print(
                "Can not locate jellyfish, please input full path of jellyfish\n"
            )

            parser.print_help()

            sys.exit(1)

    else:

        jellyfishpath = which('jellyfish')

        if jellyfishpath:

            jellyfishversion = jellyfish.jfversion(jellyfishpath[0])

            if jellyfishversion == 'None':

                print(
                    "Can not locate jellyfish, please input full path of jellyfish\n"
                )

                parser.print_help()

                sys.exit(1)

            else:

                args.jellyfish = jellyfishpath[0]

        else:

            print(
                "Can not locate jellyfish, please input full path of jellyfish\n"
            )

            parser.print_help()

            sys.exit(1)
    # End check jellyfish

    if not os.path.exists(args.genome):

        print("Can not locate close related genome file, please input a close related genome file.\n")

        parser.print_help()

        sys.exit(1)

    if not os.path.exists(args.target):

        print("Can not locate target genome file, please input target genome file.\n")

        parser.print_help()

        sys.exit(1)

    if not os.path.exists(args.input):

        print("Can not locate input file, please input input file.\n")

        parser.print_help()

        sys.exit(1)


    # Start check saved folder
    if os.path.exists(args.saved):

        print(
            args.saved,
            "exists. Everything in this folder will be remove. Press Y or N to continue: "
        )

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
    print("#" * 40)

    print("bwa version:", args.bwa, bwaversion)

    print("jellyfish version:", args.jellyfish, jellyfishversion)

    print("close related genome file:", args.genome)

    print("target genome file:", args.target)

    print("input file:", args.input)

    print("result output folder:", os.path.realpath(args.saved))

    print("threads number:", args.threads)

    print("#" * 40)

    return args


def getch():
    """
    For yes/no choice
    """
    import sys
    import tty
    import termios
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

    parser = argparse.ArgumentParser(
        description=
        "Find probes which can hybridize to a close related species.",
        prog='ChorusHomo',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n"
        "  ChorusHomo -i probe.bed -g close_related_genome.fasta -tg target_genome.fasta \\ \n"
        "             -j /opt/software/jellyfish/bin/jellyfish -b /opt/software/bwa/bwa \\ \n"
        "             -t 4 -s sample"
    )

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument(
        '-j',
        '--jellyfish',
        dest='jellyfish',
        help='The path where Jellyfish software installed')

    parser.add_argument(
        '-b',
        '--bwa',
        dest='bwa',
        help='The path where BWA software installed')

    parser.add_argument(
        '-g',
        '--genome',
        dest='genome',
        help=
        'Fasta format genome file which the probe will align to, should include all sequences from genome',
        required=True)

    parser.add_argument(
        '-tg',
        '--target',
        dest='target',
        help=
        'Fasta format genome file which the probe were generated from, should include all sequences from genome',
        required=True)

    parser.add_argument(
        '-i',
        '--input',
        dest='input',
        help=
        'BED format input file, contains oligo probes generated from Chorus2',
        required=True)

    parser.add_argument(
        '-s',
        '--save',
        dest='saved',
        help='The output folder for saving results',
        default='probes')

    parser.add_argument(
        '-t',
        '--threads',
        dest='threads',
        help='Number of threads or CPUs to use. (Default: 1)',
        default=1,
        type=int)

    return parser


if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)
