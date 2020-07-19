import pandas as pd


class Probes:

    def __init__(self, filename):

        self.filename = filename

        self.__readprobe()

    def readprobe(self):

        filename = self.filename

        probe = pd.read_table(filename, header=None)

        probe[3] = (probe[1]/1000).astype('int')

        print(probe)

        maxlength=probe[2][probe[2].idxmax(axis=1)]

        chrs = probe[0].unique().tolist()

        chrlens = dict()

        for chrnow in chrs:

            chrlen = max(probe[probe[0]==chrnow][3])

            print(chrnow, chrlen)

            chrlens[chrnow] = chrlen + 1

        self.chrs = chrs

        self.maxlength = maxlength + 1

        self.probe = probe

        self.chrlens = chrlens

    __readprobe = readprobe