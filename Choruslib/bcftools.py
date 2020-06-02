from subprocess import Popen, PIPE
from Choruslib import subprocesspath
from Choruslib import revcom
import re


def bamtobcf(bcfbin, reffile, bamfile, outbcf):
    bcfbin = subprocesspath.subprocesspath(bcfbin)

    reffile = subprocesspath.subprocesspath(reffile)

    bamfile = subprocesspath.subprocesspath(bamfile)

    outbcf = subprocesspath.subprocesspath(outbcf)

    bcfcmd = ' '.join(
        ['bcftools mpileup -E -d 500 -L 500 -Ou -f', reffile, bamfile, '|  bcftools call -cv -Ob -o', outbcf])

    print(bcfcmd)

    bcfrun = Popen(bcfcmd, shell=True)

    bcfrun.communicate()

    bcfidxcmd = ' '.join(['bcftools index', bamfile])

    print(bcfidxcmd)

    bcfidxrun = Popen(bcfidxcmd, shell=True)

    bcfidxrun.communicate()

    return True


def getconsensus(bcftoolspath, bcffile, chrom, start, end, seq, strand='+'):
    """
    get consensus by using bcftools 
    """
    bcftoolspath = subprocesspath.subprocesspath(bcftoolspath)
    bcffile = subprocesspath.subprocesspath(bcffile)
    seqlen = str(len(seq))
    pat = re.compile('[ATCG]{' + seqlen + '}')
    if strand == '-':
        seq = revcom.revcom(seq)
    fastring = '\'>' + chrom + ':' + start + '-' + end + '\\n' + seq + '\''
    bcfcon_command = ' '.join(['echo', fastring, '|' + bcftoolspath + ' consensus ', bcffile])

    p = Popen(bcfcon_command, shell=True, stdin=PIPE, stdout=PIPE)

    consensus = ''
    for i in p.stdout:
        i = i.decode('utf-8').rstrip('\n')
        #         print(i)
        if re.search(pat, i):
            consensus = i
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