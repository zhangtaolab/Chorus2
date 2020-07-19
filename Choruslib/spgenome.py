from pyfasta import Fasta
from os import path

def spgenome(fafile, outdir, maxsize=1000000000):


    spfiles = list()
    if path.exists(fafile):

        outfiles = dict()

        subfiles = dict()

        infa = Fasta(fafile)

        # nowsub = 0

        nowlen = 0

        for chrom in infa.keys():

            chrlen = len(infa[chrom])

            nowlen = nowlen+chrlen

            nowsub = int(nowlen/maxsize)

            if nowsub not in subfiles:

                subfilename = 'tmpfile' + str(nowsub) + '.fa'

                subfile = path.join(outdir,subfilename)

                spfiles.append(subfile)

                subfiles[nowsub] = open(subfile,'w')

            # outfiles[chrom] = nowsub

            print('>', chrom, sep='', file=subfiles[nowsub])

            print(infa[chrom], file=subfiles[nowsub])

        for nowsub in subfiles:

            subfiles[nowsub].close()


    else:
        print("Can't find ", fafile)

    return spfiles

if __name__ == '__main__':

    pass