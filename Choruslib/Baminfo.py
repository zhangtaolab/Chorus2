from Choruslib import openBam


class Baminfo:

    def __init__(self, bamfile):
        self.bamfile = bamfile

        self.samfile = openBam.openBam(self.bamfile)

        self.chrlen = self.getchrlen()

    def getchrlen(self):
        ref_lengths = self.samfile.lengths

        sam_ref = self.samfile.references

        refere_ncenumber = self.samfile.nreferences

        chrlen = dict()

        for i in range(refere_ncenumber):
            chr_length = ref_lengths[i]

            chrlen[sam_ref[i]] = chr_length

        return chrlen

# if __name__ == '__main__':
#
#     baminfo = Bamfile('../tests/inputfile/npstest.bam')
#
#     print()
