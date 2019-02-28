from pyfasta import Fasta


def writebed(probelist, outbedfile):

    '''probe list format:

        chr\tstart\tend

    '''
    outio = open(outbedfile, 'w')

    for pbnow in probelist:

        print(pbnow, file=outio)

    outio.close()


def writefa(genomefile, bedfile, outfile):

    fastafile = Fasta(genomefile)

    bedio = open(bedfile, 'r')

    outio = open(outfile, 'w')

    for lin in bedio.readlines():

        lin = lin.rstrip()

        chrnow, start, end = lin.split('\t')

        seqid = '>' + chrnow + ':' + start + '-' + end

        nowseq = fastafile[chrnow][int(start):int(end)]

        print(seqid, file=outio)

        print(nowseq, file=outio)

    bedio.close()

    outio.close()

    # return True



