from multiprocessing import Pool
import threading
import sys
import time
import os
import Choruslib.prefilter as prefilter
import Choruslib.primer3_filter as primer3_filter
import Choruslib.loadfa as loadfa
import Choruslib.blat as blat
import random
import Choruslib.exonerate as exonerate


def blat_getseqlen(par):

    sequence = par['sequence']

    blatport = par['blatport']

    gfcpath = par['gfcpath']

    pblength = par['pblength']

    file2bit = par['file2bit']

    dtm = par['dtm']

    minIdentity = par['minIdentity']

    st = 0

    # pblength = 48

    seqlen = len(sequence) - pblength

    print(seqlen, blatport, gfcpath, sequence[0:10])

    probelist = list()

    while st < seqlen:

            sp = st+pblength

            nowseq = sequence[st:sp]

            st += 5

            if not nowseq:

                continue

            if not prefilter.atcg_filter(nowseq, 7):

                if 0 < blat.blat_search_sequence(gfcpath, sequence=nowseq, blatport=blatport, minIdentity=minIdentity, file2bit=file2bit) < 2:

                    if not primer3_filter.primer3_filter(sequence=nowseq, dtm=dtm):

                        probelist.append(nowseq)

                        st += pblength

    return probelist


def blat_pbsearch(fafile, file2bit, gfspath, gfcpath, threadnum, pblength, setpsize, dtm=10, minIdentity=75):

    fa = loadfa.loadfa(fafile)

    fa = fa.replace('n', '')

    fa = fa.replace('N', '')

    windlen = 1000000 # for real
    # windlen = 10000 #for test
    seqlen = len(fa)

    sptime = int(seqlen/windlen)

    pars = list()

    pt = 0

    chp = int(threadnum/6 - 1/6) + 1

    for por in range(chp):

        blatport = por+10010

        t = threading.Thread(target=blat.start_gfServer, kwargs={'file2bit':file2bit, 'gfspath':gfspath, 'stepsize':setpsize, 'blatport':blatport})

        t.daemon = True

        t.start()

    print("starting gfServer ... ")

    time.sleep(100)

    for sc in range(sptime+1):

        nowstart = sc * windlen

        nowend = nowstart + windlen + 50

        if nowend > seqlen:

            nowend = seqlen-1

        nowseq = fa[nowstart:nowend]



        blatport = 10010 + random.randint(0, chp-1)

        pt = 1 + pt

        par = dict()

        print (nowstart,nowend,blatport)

        par['sequence'] = nowseq

        par['blatport'] = blatport

        par['gfcpath'] = gfcpath

        par['pblength'] = pblength

        par['file2bit'] = file2bit

        par['dtm'] = dtm

        par['minIdentity'] = minIdentity

        pars.append(par)

    par2 = dict()

    blatport = 10010

    nowstart = windlen * sptime

    nowend = seqlen-1

    nowseq = fa[nowstart:nowend]

    par2['sequence'] = nowseq

    par2['blatport'] = blatport

    par2['gfcpath'] = gfcpath

    par2['pblength'] = pblength

    par2['file2bit'] = file2bit

    par2['dtm'] = dtm

    par2['minIdentity'] = minIdentity

    pars.append(par2)

    pool = Pool(threadnum)

    res = pool.map(blat_getseqlen, pars)

    probe = dict()

    for i in res:

        for j in i:

            if j in probe:

                probe[j] = +1
            else:
                probe[j] = 1


    blat.stop_gfServer()

    # probecount = Counter
    #
    # for probeseq in probe:
    #
    #     probecount[probeseq] += 1

    probelist = list()

    for probeseq in probe:

        if probe[probeseq] == 1:

            probelist.append(probeseq)

    return probelist


def exonerate_pbsearch(fafile, esifile, espath, ep, threadnum, pblength, setpsize, dtm=10, minIdentity=75,
                       exonerateport=10005):

    fa = loadfa.loadfa(fafile)

    fa = fa.replace('n', '')

    fa = fa.replace('N', '')

    seqlen = len(fa)

    windlen = 10000

    sptime = int(seqlen/windlen)+1

    pars = list()

    #start exonerate server

    # exonerate.start_exonerate_server(espath=espath, esifile=esifile, threadnum=24,  exonerateport=exonerateport)

    time.sleep(10)

    for sc in range(sptime):

        nowstart = sc * windlen

        nowend = nowstart + windlen + 50

        if nowend > seqlen:

            nowend = seqlen-1

        nowseq = fa[nowstart:nowend]

        # print(nowseq)

        par = dict()

        par['sequence'] = nowseq

        par['exonerateport'] = exonerateport

        par['ep'] = ep

        par['minIdentity'] = minIdentity

        par['pblength'] = pblength

        pars.append(par)

    pool = Pool(threadnum)

    res = pool.map(exonerate_getseqlen, pars)

    probe = dict()

    for i in res:

        for j in i:

            if j in probe:

                probe[j] = +1
            else:
                probe[j] = 1


    exonerate.stop_exonerate_server()


    probelist = list()

    print(probe)

    for probeseq in probe:

        if probe[probeseq] == 1:

            probelist.append(probeseq)

    return probelist


def exonerate_getseqlen(par):

    sequence = par['sequence']

    exonerateport = par['exonerateport']

    ep = par['ep']

    minIdentity = par['minIdentity']

    pblength = par['pblength']

    print(sequence[1:10], exonerateport, ep, minIdentity, pblength)

    st = 0

    # pblength = 48

    seqlen = len(sequence) - pblength

    probelist = list()

    tmpfa = sequence[1:10] + str(random.randint(1, 100)) + str(int(time.time())) + '.fa'

    tmpfaio = open(tmpfa,'w')

    tmpseq = dict()

    while st < seqlen:

            sp = st+pblength

            nowseq = sequence[st:sp]

            st += 5

            # print(nowseq)

            if not nowseq:

                continue

            if not prefilter.atcg_filter(nowseq, 7):

                if nowseq in tmpseq:

                    tmpseq[nowseq] += 1

                else:

                    tmpseq[nowseq] = 1

            for nowseq in tmpseq:

                if tmpseq[nowseq] == 1:

                    print('>',nowseq, sep='', file=tmpfaio)

                    print(nowseq, file=tmpfaio)

    tmpfaio.close()

    exonerate.exonerate_search_sequence2(ep=ep, sequencefile=tmpfa, exonerateport=exonerateport,
                                                             minIdentity=minIdentity)
    return probelist

if __name__ == "__main__":

    try:
        prb = exonerate_pbsearch(fafile='/Users/Forrest/Documents/Project/Chorus/Test/PGSC_DM_v4.03_pseudomolecules_genes_exon_test.fa',
                       esifile='/Users/Forrest/Documents/Project/Chorus/Test/Gy14.sei',
                       espath='/opt/local/bin/exonerate-server', ep='/opt/local/bin/exonerate',
                       setpsize=5,
                       exonerateport=10005,
                       threadnum=2,
                       pblength = 45, dtm=10, minIdentity=75)

        print(prb)
    except KeyboardInterrupt:
        exonerate.stop_exonerate_server()
        sys.stderr.write("User interrupt\n")

        sys.exit(0)