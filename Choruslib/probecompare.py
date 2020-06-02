def compareseq(seq1, seq2):
    len1 = len(seq1)
    len2 = len(seq2)
    dif = '-1'
    same = 0
    if len1 == len2:

        for i in range(len1):

            if seq1[i] == seq2[i]:
                same += 1

        dif = len1 - same

    else:
        dif = '-1'

    return dif


def getconsensusprobe(probelist):

    probelen0 = len(probelist[0])

    consensusprobe = ''

    for probenow in probelist:

        if probelen0 != len(probenow):
            consensusprobe = -1
            consensusscore = -1
            break

    if not consensusprobe == -1:

        consensuslist = list()
        consensuscount = list()
        for i in range(len(probelist[0])):

            nowdict = dict()

            for j in range(len(probelist)):
                #                 print(i,j)
                if probelist[j][i] in nowdict:

                    nowdict[probelist[j][i]] += 1

                else:

                    nowdict[probelist[j][i]] = 1

            consensuslist.append(max(nowdict.items(), key=operator.itemgetter(1))[0])

            consensuscount.append(str(max(nowdict.items(), key=operator.itemgetter(1))[1]))

        consensusprobe = ''.join(consensuslist)
        consensusscore = ''.join(consensuscount)

    #     print("indef", consensusprobe)

    return (consensusprobe, consensusscore)
