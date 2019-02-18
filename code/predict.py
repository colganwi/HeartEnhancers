import numpy as np
from Bio import SeqIO
from file_loading import *
import matplotlib.pyplot as plt


alignments = load_alignments("data/alignments.txt")

for alignment in alignments[0:10]:
    print(alignment[0]+'\n')

    cr = alignment[1]
    cs = alignment[2]
    SIZE = len(cr)
    cr_ggaw = np.zeros(SIZE,np.int8)
    cs_ggaw = np.zeros(SIZE,np.int8)
    cr_atta = np.zeros(SIZE,np.int8)
    cs_atta = np.zeros(SIZE,np.int8)
    cr_tgtt = np.zeros(SIZE,np.int8)
    cs_tgtt = np.zeros(SIZE,np.int8)
    cr_coding = np.zeros(SIZE,np.int8)
    cs_coding = np.zeros(SIZE,np.int8)
    identity = np.zeros(SIZE,np.float16)
    ggaw_amount = np.zeros(SIZE,np.float16)
    atta_amount = np.zeros(SIZE,np.float16)
    tgtt_amount = np.zeros(SIZE,np.float16)
    for i in range(SIZE-4):
        cr_mer = cr[i:i+4]
        cs_mer = cs[i:i+4]
        if cr_mer == 'GGAT' or cr_mer == 'GGAA' or cr_mer == 'ATCC' or cr_mer == 'TTCC':
            cr_ggaw[i] = 1
        if cs_mer == 'GGAT' or cs_mer == 'GGAA' or cs_mer == 'ATCC' or cs_mer == 'TTCC':
            cs_ggaw[i] = 1
        if cr_mer == 'ATTA' or cr_mer == 'TAAT':
            cr_atta[i] = 1
        if cs_mer == 'ATTA' or cs_mer == 'TAAT':
            cs_atta[i] = 1
        if cr_mer == 'TGTT' or cr_mer == 'AACA':
            cr_tgtt[i] = 1
        if cs_mer == 'TGTT' or cs_mer == 'AACA':
            cs_tgtt[i] = 1

    conserved_ggaw = cr_ggaw * cs_ggaw
    conserved_atta = cr_atta * cs_atta
    conserved_tgtt = cr_tgtt * cs_tgtt

    for i in range(SIZE):
        if cr[i] == '#':
            cr_coding[i] = 1

    for i in range(SIZE):
        if cs[i] == '#':
            cs_coding[i] = 1

    region = []
    for i in range(150):
        if cr[i] == cs [i] and cr[i] != '-' and cr[i] != '#':
            region += [1]
        else:
            region += [0]
    for i in range(150,SIZE):
        del region[0]
        if cr[i] == cs [i] and cr[i] != '-' and cr[i] != '#':
            region += [1]
        else:
            region += [0]

        score = sum(region)/25.0 + sum(cr_ggaw[i-150:i]) + 2*sum(cs_ggaw[i-150:i]) + 2*sum(cr_atta[i-150:i]) + sum(cs_atta[i-150:i]) + sum(cr_tgtt[i-150:i]) + sum(cs_tgtt[i-150:i])
        identity[i] = sum(region)/150.0
        ggaw_amount[i] = sum(conserved_ggaw[i-150:i])/10.0
        atta_amount[i] = sum(conserved_atta[i-150:i])/10.0
        tgtt_amount[i] = sum(conserved_tgtt[i-150:i])/10.0

    score = ((identity-.5) + (ggaw_amount-.1) + (atta_amount-.1) + (tgtt_amount-.1))
    score[score < 0] = 0
    #score[10000:] = 0
    if(max(score > .2)):
        i = np.argmax(score)
        print(str(i-150)+": "+cr[i-150:i])
        print(str(i-150)+": "+cs[i-150:i])

    fig, ax = plt.subplots()
    ax.plot(identity,label='identity')
    ax.plot(ggaw_amount,label='ggaw')
    ax.plot(atta_amount,label='atta')
    ax.plot(tgtt_amount,label='tgtt')
    ax.plot(score,color='black',label='score')
    ax.fill_between(range(SIZE),cr_coding,color='r',label='coding',alpha=0.5,linestyle='None')
    ax.fill_between(range(SIZE),cs_coding,color='b',label='coding',alpha=0.5,linestyle='None')
    ax.legend(loc='upper right')
    plt.show()
