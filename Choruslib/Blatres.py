class Blatres:

    def __init__(self, seq, blatlines):

        self.seq = seq

        self.blatlines = blatlines

        self.matchednumber = 0

        self.identicnumber = 0

        self.blatidenticlist = list()

        self._processblatline()

    def processblatline(self):

        for lin in self.blatlines:

            lin = lin.rstrip("\n")

            (match,mismatch,repmatch,Ns, qgapcount, qgapbase, tgapcount, tgapbases,
             strand, qname, qsize, qstart, qend, tname, tsize, tstart,tend,blockcount,
             blocksize,qStarts,qStarts) = lin.split("\t")

            nscore = int(mismatch)+int(repmatch)+int(qgapbase)+int(tgapbases)

            self.matchednumber += 1

            if nscore == 0:

                blatidenticres = dict()

                blatidenticres[qname] = qname

                blatidenticres[tname] = tname

                blatidenticres[tstart] = tstart

                blatidenticres[tend] = tend

                self.blatidenticlist.append(blatidenticres)

                self.identicnumber += 1

    _processblatline = processblatline