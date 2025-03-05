import os
import sys


class FastaNotFound(Exception):

    pass


def loadfa(fastafile):

    if not os.path.exists(fastafile):
            raise FastaNotFound('"' + fastafile + '"')

    fh = open(fastafile, 'r')

    seqs = ''

    for line in fh:

        line = line.rstrip()
        # print (line)
        if not line:

            continue

        if line[0] == ">":

            pass

        else:

            seqs = seqs + line

    fh.close()

    return seqs


# if __name__ == "__main__":
#
#     try:
#
#         fa = loadfa(fastafile='/Users/Forrest/Documents/Project/Chorus/Test/PGSC_DM_v4.03_pseudomolecules_genes_exon_test.fa')
#
#         print(fa)
#
#     except KeyboardInterrupt:
#
#         sys.stderr.write("User interrupt\n")
#
#         sys.exit(0)