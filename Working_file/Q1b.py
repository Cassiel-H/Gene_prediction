# Jiachen Huo
# 260779330

from importlib import import_module
import sys
from Bio import SeqIO
import numpy as np
import math as mt


def main():

    # Note all following assume state [Inter, Start, Middle, Stop] is in order of [0,1,2,3]

    # .fa file and configuration file
    seqfn = sys.argv[1]
    config = sys.argv[2]

    # load the configuration file as import_module
    temp = import_module(config.split('.')[0])
    inputFig = temp.config()

    # get merged nucleotide sequence as a string
    merged = merge(seqfn)

    # write out the file by iterating through all .fa file's record
    records = SeqIO.parse(open(seqfn), 'fasta')
    output = open(r"1b_out.gff3", "w")
    for i in records:
        grid, track, l, seq_id= initialize(inputFig, i)  # implement viterbi as HMM
        #l, seq_id = traceBack(grid, track, i)  # trace-back procedure
        rec = []
        for j in range(len(l) - 1, -1, -1):
            rec += [[seq_id, 'ena', 'CDS', str(l[j][0]), str(l[j][1])]]
        for x in rec:
            output.writelines(str('   '.join(x)) + '\n')


# apply logarithm to both emission and transition probability to avoid underflow
# math domain error could happen when apply log to base 0 or near 04
def log(x):
    return mt.log(x) if x > 0 else -mt.inf


# merge the given sequence as a concatenated string
def merge(seqfn):
    sequence = SeqIO.parse(open(seqfn), 'fasta')
    genome = ''.join(str(i.seq) for i in sequence)
    return genome


# viterbi algorithm of HMM
# take input as 1.the config class object which contains e_table and t_table; 2. the
def initialize(inputFig, sequence):
    seqfn = str(sequence.seq)
    length = len(seqfn)
    grid = np.full((4, length), -mt.inf)  # create nd-array for dynamic programming table, initial set to -inf
    track = np.full((4, length), -1)  # create same-sized nd-array for a trace back pointer, initial -1

    grid[0, 0] = mt.log(inputFig.initial[0]) + mt.log(inputFig.emission[0][seqfn[0]])  # initialize the grid
    grid[1, 0] = grid[2, 0] = grid[3, 0] = 0
    size = [1, 3, 3, 3]

    # The dynamic procedure
    for j in range(length):    # loop all column (with length=input sequence)
        for i in range(4):      # loop all four state
            store = -mt.inf
            pointer = -1
            x = seqfn[j - size[i] + 1:j + 1]
            for p in range(4):  # as j-1 column's states
                check = max(j - size[p], 0)
                temp = grid[p, check] + log(inputFig.transition[p][i])  # temporary probability if previous state was p's position
                if temp > store:
                    store = temp
                    pointer = p

            grid[i, j] = store + log(inputFig.emission[i][x])
            track[i, j] = pointer

    # Trace back procedure
    gene_region = []
    index = int(np.argmax(grid[:, -1]))                 # index refers to the index of the row of previous column,
    region = False                                      # that gives the largest value to this column
    start = end = 0
    for i in range(len(track[1]) - 1, -1, -1):         # loop from last column of track
        if index == 3 and region is False:             # if we obtain index=3 (stop state), and indicator region= False
            end = i + 1                                # means we will enter the predicted gene, set end point= i+1
            region = True
        if index == 1 and region is True:              # if obtain index=1 (Start state), and region = True
            start = i - 3                              # Means we reach the start of one gene, set start point= i-3
            region = False
            gene_region += [[start, end]]
        index = track[index, i]

    return grid, track, gene_region, sequence.id       # return the currently reversed gene_region list of list,
                                                       # and seq.id


if __name__ == '__main__':
    main()
