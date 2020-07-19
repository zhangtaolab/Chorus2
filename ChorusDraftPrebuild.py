from pyfasta import Fasta
import sys
import argparse
import os


def check_options(parser):

    args = parser.parse_args()

    if args.input:

        if not os.path.exists(args.input):

            print("Can not locate input file, please input input file.\n")

            parser.print_help()

            sys.exit(1)

    return args


def get_options():

    parser = argparse.ArgumentParser(description="Combine short sequence to speed up oligo search", prog="ChorusDraftPrebuild")

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('-i', '--input', dest='input', help='Fasta format input file contains short sequences',
                        required=True, type=str)

    parser.add_argument('-o', '--output', dest='output', help='Fasta format output file with combined long sequences for speeding up oligo search. (default: output.fa)',
                        default='output.fa', type=str)

    return parser


def main():

    args = check_options(get_options())

    fain = Fasta(args.input)

    faout = open(args.output, 'w')

    minlen = int(1e6)

    print(minlen)

    shortseq = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'

    breacker = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'

    shortlist = list()

    for chrome in fain.keys():

        if len(fain[chrome]) < minlen:
            # print(chrome, len(fain[chrome]))
            # shortseq = shortseq + str(fain[chrome]) + breacker
            shortlist.append(chrome)

        else:
            print(chrome, len(fain[chrome]))
            print('>%s' % chrome, file=faout)
            print(fain[chrome], file=faout)

    print('>shortsequences', file=faout)

    for chrome in shortlist:

        print(str(fain[chrome]),shortseq,sep='',end='', file=faout)

#    print(shortseq, file=faout)


    faout.close()


if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)
