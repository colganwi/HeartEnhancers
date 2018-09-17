"""
These functions are for file loading
"""
import numpy as np
from Bio import SeqIO
import gffutils

def load_list(filepath):
    file= open(filepath,"r")
    lines = file.readlines()
    list = []
    for line in lines:
        list.append(line.strip('\n'))
    file.close
    return list

def load_genes(filepath):
    file= open(filepath,"r")
    lines = file.readlines()
    genes = []
    for line in lines:
        genes.append(line.strip('\n').split('\t'))
    file.close
    return genes

def load_blast(filepath):
    file= open(filepath,"r")
    lines = file.readlines()
    hits = {}
    for line in lines:
        line = line.split()
        hits[line[0]] = line[1]
    file.close
    return hits

def load_gff(filepath):
    return gffutils.FeatureDB(filepath, keep_order=True)

def load_fasta(filepath):
    reads = SeqIO.parse(open(filepath),'fasta')
    return reads

def load_alignments(filepath):
    file= open(filepath,"r")
    lines = file.readlines()
    alignments = []
    alignment = ['','']
    gene = ''
    first = True
    i = 0
    while i < len(lines):
        if lines[i][:2] == 'KH':
            if not first:
                alignments += [[gene,alignment[0],alignment[1]]]
            first = False
            gene = lines[i].split(' ')[0]
            alignment = ['','']
            i+=2
        else:
            alignment[0] += lines[i].strip('\n')
            alignment[1] += lines[i+1].strip('\n')
            i+=3
    return alignments
