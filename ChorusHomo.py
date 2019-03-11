import os
import sys
from Choruslib import bwa
import Chorus


def main():
    # get args
    args = Chorus.check_options(Chorus.get_options())
    print(args)

    # init input and output files
    related_genome = args.genome
    probe_file = args.input
    probe_fasta = bed_to_fa(probe_file, args.saved)
    homo_probe = make_saved_file(args.saved, probe_file, "homo.bed")

    sys.stderr.write("input probe {}\n".format(probe_file))

    # align probes to genome
    if not check_bwa_index(related_genome, args.bwa, args.saved):
        sys.stderr.write("Failed to get index for ", related_genome, "\n")
        sys.exit(0)

    bwafiltedpb = bwa.bwafilter(
        bwabin=args.bwa,
        reffile=related_genome,
        inputfile=probe_fasta,
        minas=args.length,
        maxxs=int(args.length * args.homology / 100),
        threadnumber=args.threads)
    #
    hom_probe_dict = dict()
    for op in bwafiltedpb:
        seq, chr, st = op.strip().split()
        head = chr + st
        hom_probe_dict[head] = 1


# add homology information after each oligos
    with open(homo_probe, "w") as o:
        with open(probe_file, "r") as i:
            for line in i:
                arr = line.strip().split()
                hh = arr[0] + arr[1]
                if hh in hom_probe_dict:
                    arr.append("H")
                else:
                    arr.append("N")
                save_probe(arr, o)

    sys.stderr.write("record homologous oligos\n")


def make_saved_file(saved_path, file, extention):
    ''' add full path to file '''
    name = os.path.basename(file)
    base = os.path.splitext(name)[0]
    saved_file = "".join((saved_path, "/", base + "." + extention))
    return saved_file


def check_bwa_index(related_genome, bwa, saved_path):
    ''' check if bwa index exist '''
    index = related_genome + ".sa"
    if os.path.isfile(index):
        sys.stderr.write("index already existed\n")
        return True
    else:
        sys.stderr.write("make index file for genome\n")
        bwa.bwaindex(bwa, related_genome, saved_path)
        return True

    return False


def bed_to_fa(probe_file, saved_path):
    ''' convert bed to fasta, header is the position of the seq.
    All other information in the bed file is lost'''
    fasta = make_saved_file(saved_path, probe_file, "fa")
    with open(fasta, "w") as o:
        with open(probe_file, "r") as i:
            for line in i:
                chr, st, ed, seq = line.strip().split()
                header = "_".join((chr, st, ed))
                o.write(">" + header + "\n" + seq + "\n")

    return fasta


def save_probe(one_oligo, output_handle):
    ''' save oligo to file '''
    output_handle.write("\t".join(one_oligo))
    output_handle.write("\n")


if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)
