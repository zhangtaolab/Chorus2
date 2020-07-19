import re


def atcg_filter(sequence, rept=7):

    nmatch = re.compile(r'[Nn]')

    nfilter = bool(nmatch.search(sequence))

    atcgmatch = re.compile(r'[aA]{%s}|[tT]{%s}|[cC]{%s}|[Gg]{%s}|[Nn]{%s}' % (rept, rept, rept, rept, rept))

    atcgfilter = bool(atcgmatch.search(sequence))

    filterseq = nfilter or atcgfilter

    return filterseq


def gc_filter(sequence, maxgc, mingc):

    pass
