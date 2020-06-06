from subprocess import Popen, PIPE
from Choruslib import subprocesspath
from Choruslib import revcom
import re


class BcfConsensusRuner():

    def __init__(self, probestr, bcftoolspath, bcffile, sample):

        self.probestr = probestr

        self.bcftoolspath = bcftoolspath

        self.bcffile = bcffile

        self.sample = sample

def bamtobcf(bcfbin, reffile, bamfile, outbcf):
    bcfbin = subprocesspath.subprocesspath(bcfbin)

    reffile = subprocesspath.subprocesspath(reffile)

    bamfile = subprocesspath.subprocesspath(bamfile)

    outbcf = subprocesspath.subprocesspath(outbcf)

    bcfcmd = ' '.join(
        [bcfbin,' mpileup -E -d 500 -L 500 -Ou -f', reffile, bamfile, '| ', bcfbin,' call -cv -Ob -o', outbcf])

    print(bcfcmd)

    bcfrun = Popen(bcfcmd, shell=True)

    bcfrun.communicate()

    bcfidxcmd = ' '.join([bcfbin, ' index', outbcf])

    print(bcfidxcmd)

    bcfidxrun = Popen(bcfidxcmd, shell=True)

    bcfidxrun.communicate()

    return True

def probestrtoconsensus(bcfconsensusruner):

    (chrom, start, end, seq, score, strand) = str(bcfconsensusruner.probestr).rstrip().split("\t")

    if strand == '-':
        seq = revcom.revcom(seq)

        strand = '+'

    consensusprobe = getconsensus(bcftoolspath=bcfconsensusruner.bcftoolspath,
                                  bcffile=bcfconsensusruner.bcffile,
                                  chrom=chrom,
                                  start=start,
                                  end=end,
                                  seq=seq,
                                  sample=bcfconsensusruner.sample
                                  )
    res = dict()

    res['probestr'] = bcfconsensusruner.probestr

    res['consensusprobe'] = consensusprobe

    return res

def getconsensus(bcftoolspath, bcffile, chrom, start, end, seq, sample, strand='+'):
    """
    get consensus by using bcftools 
    """
    bcftoolspath = subprocesspath.subprocesspath(bcftoolspath)
    bcffile = subprocesspath.subprocesspath(bcffile)
    seqlen = str(len(seq))
    pat = re.compile('[ATCG]{' + seqlen + ',}')
    if strand == '-':
        seq = revcom.revcom(seq)
    fastring = '\'>' + chrom + ':' + start + '-' + end + '\\n' + seq + '\''
    bcfcon_command = ' '.join(['echo', fastring, '|' + bcftoolspath + ' consensus -s', sample, bcffile])

    consensus = 'N'*len(seq)

    try:
        p = Popen(bcfcon_command, shell=True, stdin=PIPE, stdout=PIPE)

        for i in p.stdout:
            i = i.decode('utf-8').rstrip('\n')
            #         print(i)
            if pat.search(i):
                consensus = pat.search(i)[0]
    except:
        print("warnning: ", bcfcon_command, " ##")
    #             print('c:',consensus)
    return str(consensus)


def bcftoolsversion(bcftoolsbin):
    """

    :param bcftoolsbin: bwa bin path
    :return: string, version of bcftools
    """

    bcftoolscmd = [bcftoolsbin]

    bcftoolsrun = Popen(bcftoolscmd, stdout=PIPE, stderr=PIPE)

    pat = re.compile('Version')

    version = 'None'

    for i in bcftoolsrun.stderr.readlines():

        i = i.decode('utf-8').rstrip('\n')

        if re.search(pat, i):

            version = i

    bcftoolsrun.communicate()

    return version