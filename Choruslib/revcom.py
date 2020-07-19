def revcom(sequence):

    revdic = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','a':'t','t':'a','c':'g','g':'c'}

    return "".join([revdic[base] for base in reversed(sequence)])