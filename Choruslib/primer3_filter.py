import primer3


# def primer3_filter(sequence, mintm=0, maxhtm=0, dtm=10):
#
#     primer3ft = True
#
#     tm = primer3.calcTm(sequence)
#
#     htm = primer3.calcHairpinTm(sequence)
#
#     if mintm*maxhtm == 0:
#
#         if (tm - htm) > dtm:
#
#             primer3ft = False
#
#     else:
#
#         pass
#
#     return primer3ft

def revcom(sequence):

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    reverse_complement = "".join(complement.get(base, base) for base in reversed(sequence))

    return reverse_complement


def primer3_filter(sequence, mintm=37, maxhtm=35, dtm=10):

    primer3ft = False

    tm = primer3.calcTm(sequence)

    htm = primer3.calcHairpinTm(sequence)

    if tm < mintm:

        primer3ft = True

    if htm > maxhtm:

        primer3ft = True

    if (tm-htm) < dtm:

        primer3ft = True

    # print(sequence, tm, htm, dtm)

    return primer3ft


def primer3_filter_withRprimer(sequence, rprimer, mintm=37, maxhtm=35, dtm=10):

    primer3ft = False

    tm = primer3.calcTm(sequence)

    fseq = rprimer + sequence

    htmF = primer3.calcHairpinTm(fseq)

    rseq = rprimer + revcom(sequence)

    htmR = primer3.calcHairpinTm(rseq)

    if tm < mintm:

        primer3ft = True

    if htmF > maxhtm:

        primer3ft = True

    if (tm-htmF) < dtm:

        primer3ft = True

    if htmR > maxhtm:

        primer3ft = True

    if (tm-htmR) < dtm:

        primer3ft = True
    # print(sequence, tm, htm, dtm)

    return primer3ft


def primer3_cal(sequence, mintm=37, maxhtm=37, dtm=10):

    primer3ft = True

    tm = primer3.calcTm(sequence)

    htm = primer3.calcHairpinTm(sequence)

    if tm < mintm:

        primer3ft = False

    if htm > maxhtm:

        primer3ft = False

    if (tm-htm) > dtm:

        primer3ft = False

    return (sequence, primer3ft)
