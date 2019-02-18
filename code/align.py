import gffutils
import numpy as np
from Bio import SeqIO
from Bio import pairwise2
from file_loading import *

SIZE = 10000
genes = load_genes("data/heart_genes.txt")
cr_genes = load_gff("data/cr_genes.db")
cs_genes = load_gff("data/cs_genes.db")
cr_genome = load_fasta("data/cr_genome.fasta")
cs_genome = load_fasta("data/cs_genome.fasta")
cr_seqs = {}
cs_seqs = {}
for seq in cr_genome:
    cr_seqs[seq.id] = seq.seq
for seq in cs_genome:
    cs_seqs[seq.id] = seq.seq

file= open("data/alignments.txt","w")
for gene in genes:

    file.write(gene[0]+' '+gene[1]+'\n\n')
    print(gene[0]+' '+gene[1]+'\n')

    cr_gene = cr_genes[gene[0]]
    cs_gene = cs_genes[gene[1]]

    length = 0
    if(cr_gene.strand == '+'):
        cr_start = cr_gene.start
    else:
        cr_start = cr_gene.end

    cr_mask = np.zeros(SIZE*2,np.int8)
    for CDS in cr_genes.region(region=(cr_gene.seqid,cr_start-SIZE,cr_start+SIZE), featuretype=['CDS']):
        cr_mask[CDS.start-cr_start-SIZE:CDS.end-cr_start-SIZE] = 1
    if(cr_gene.strand == '+'):
        cr_upstream = cr_seqs[cr_gene.chrom][cr_start-SIZE:cr_start]
        cr_downstream = cr_seqs[cr_gene.chrom][cr_start:cr_start+SIZE]
    else:
        cr_upstream = cr_seqs[cr_gene.chrom][cr_start:cr_start+SIZE].reverse_complement()
        cr_downstream = cr_seqs[cr_gene.chrom][cr_start-SIZE:cr_start].reverse_complement()
        cr_mask = np.flip(cr_mask,0)

    clength = 0
    if(cs_gene.strand == '+'):
        cs_start = cs_gene.start
    else:
        cs_start = cs_gene.end

    cs_mask = np.zeros(SIZE*2,np.int8)
    for CDS in cs_genes.region(region=(cs_gene.seqid,cs_start-SIZE,cs_start+SIZE), featuretype=['exon']):
        cs_mask[CDS.start-cs_start-SIZE:CDS.end-cs_start-SIZE] = 1
    if(cs_gene.strand == '+'):
        cs_upstream = cs_seqs[cs_gene.chrom][cs_start-SIZE:cs_start]
        cs_downstream = cs_seqs[cs_gene.chrom][cs_start:cs_start+SIZE]
    else:
        cs_upstream = cs_seqs[cs_gene.chrom][cs_start:cs_start+SIZE].reverse_complement()
        cs_downstream = cs_seqs[cs_gene.chrom][cs_start-SIZE:cs_start].reverse_complement()
        cs_mask = np.flip(cs_mask,0)

    alignments_upstream = pairwise2.align.globalms(cr_upstream, cs_upstream, 1, -1, -10, -.5)
    alignments_downstream = pairwise2.align.globalms(cr_downstream, cs_downstream, 1, -1, -10, -.5)

    if len(alignments_upstream) >0 and len(alignments_downstream) >0 :
        cr_alignment = alignments_upstream[0][0] + alignments_downstream[0][0]
        cs_alignment = alignments_upstream[0][1] + alignments_downstream[0][1]

        new_cr_alignment = ''
        j = 0
        for i in range(len(cr_alignment)):
            if cr_alignment[i] != '-':
                if cr_mask[j] == 1:
                    new_cr_alignment += '#'
                else:
                    new_cr_alignment += cr_alignment[i]
                j += 1
            else:
                new_cr_alignment += cr_alignment[i]
        cr_alignment = new_cr_alignment

        new_cs_alignment = ''
        j = 0
        for i in range(len(cs_alignment)):
            if cs_alignment[i] != '-':
                if cs_mask[j] == 1:
                    new_cs_alignment += '#'
                else:
                    new_cs_alignment += cs_alignment[i]
                j += 1
            else:
                new_cs_alignment += cs_alignment[i]
        cs_alignment = new_cs_alignment

    cr_alignment = str(cr_upstream) + str(cr_downstream)
    for i in range(0,2*SIZE-1,80):
        file.write(cr_alignment[i:i+80]+'\n')
        file.write(cs_alignment[i:i+80]+'\n\n')
file.close()
