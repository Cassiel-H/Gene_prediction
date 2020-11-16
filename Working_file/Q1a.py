# Jiachen Huo
# 260779330

import gffutils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pip._vendor.msgpack.fallback import xrange
import pandas as pd


class Q1a:
    def __init__(self):
        fn = 'Vibrio_cholerae.GFC_11.37.gff3'
        seqfn = "Vibrio_cholerae.GFC_11.dna.toplevel.fa"
        dbfn = "Vibrio_cholerae.GFC_11.37.db"

        db = gffutils.create_db(fn, dbfn, force=True, keep_order=True, merge_strategy='merge',
                                sort_attribute_values=True)
        db = gffutils.FeatureDB(dbfn, keep_order=True)

        self.totGenicLen, self.genicRegion, self.avgenicLen = genicSeq_avgLength(seqfn, db)
        self.totInterLen, self.intergenicRegion, self.avgInterLen = interSeq_avgLength(self.genicRegion,
                                                                                       self.totGenicLen)
        self.genic_sequence = Sequence(seqfn, self.genicRegion, "Gene")
        self.intergenic_sequence = Sequence(seqfn, self.intergenicRegion, "intergenic")

        self.list_genic_seq = list(self.genic_sequence)
        #print(self.list_genic_seq)
        self.list_inter_seq = list(self.intergenic_sequence)

        #print(list(self.list_genic_seq))
        self.merged_inter = merge(self.list_inter_seq)
        self.merged_genic = merge(self.list_genic_seq)
        #print(self.merged_genic)

        #print(self.list_genic_seq[0]+"\n"+"\n")

        # frequency table variables
        self.nuc_frequency = nuc_frenq(self.merged_inter, self.totInterLen)
        self.totalCodon = total_cod(self.merged_genic, self.totGenicLen)
        self.startCodon, self.startTable = start_cod(self.list_genic_seq, self.merged_genic)
        self.middleCodon = middle_cod(self.merged_genic, self.totGenicLen)
        self.stopCodon, self.stopTable = stop_cod(self.list_genic_seq, self.merged_genic)

        # generating sequence file as gene1-gene1887 / intergenic1-intergenic1887
        file_genic(self.list_genic_seq)
        file_intergenic(self.list_inter_seq)


# return a list of the index of genic region [[xx,xx], [xx,xx], [xx,xx] .... ];
# calculate average length of genic region in print statement;
def genicSeq_avgLength(seqfn, db):
    sequence = SeqIO.parse(open(seqfn), 'fasta')
    genic_regions = []
    counter = 0
    sumlength = 0
    for i in sequence:
        for fig in db.region(seqid=i.name, featuretype='CDS', strand='+'):
            genic_regions += [[counter + fig.start, counter + fig.end]]
            sumlength += fig.end - fig.start
        counter += len(i)
    avg_genic_length = round(sumlength / len(genic_regions))
    #print(genic_regions)
    return sumlength, genic_regions, avg_genic_length


# return a list of the index of intergenic region [[xx,xx], [xx,xx], [xx,xx] .... ];
# calculate average length of intergenic region in print statement;
def interSeq_avgLength(genicRegion, total_genic_length):
    total_intergenic_length = 0
    intergenic_seq = []
    count = 0
    if genicRegion[0][0] > 1:
        intergenic_seq += [[1, genicRegion[0][0] - 1]]
        total_intergenic_length += genicRegion[0][0] - 1
        count += 1
    for i in range(len(genicRegion) - 1):
        intergenic_seq += [[genicRegion[i][1] + 1, genicRegion[i + 1][0] - 1]]
        total_intergenic_length += genicRegion[i + 1][0] - 1 - genicRegion[i][1]
        count += 1
        if genicRegion[-1][1] < total_genic_length:
            intergenic_seq += [[genicRegion[-1][1] + 1, total_genic_length]]
            total_intergenic_length += total_genic_length - genicRegion[-1][1]
            count += 1
    avg = round(total_intergenic_length / count)
    return total_intergenic_length, intergenic_seq, avg


# return a generator of actual SeqRecord of sequence of genic/intergenic region;
# [SeqRecord(seq=seq('ATCGTGCATG.....TTT'), id=xxx, name=xxx,..
#                description='length=xxx, start:xxx, end:xxx'),
#  SeqRecord(seq=seq('TGCATGGGGC.....TTT'), id=xxx, name=xxx,..
#                description='length=xxx, start:xxx, end:xxx'),
#                               ....... ]
def Sequence(seqfn, sequenceRegion, type):
    records = SeqIO.parse(open(seqfn), 'fasta')
    seq = ''.join(str(record.seq) for record in records)

    #print(seq)
    """
    all = list(seq)
    count = 0
    for i in range(0,653):
        print(all[i])
        count += 1
    print(count)
    """

    # get sequence for genic region without concatenation
    if type == "Gene":
        sequence = (SeqRecord(
            Seq(seq[gene[0] - 1:gene[1]]),
            id="gene " + str(i + 1), description=f'length= {gene[1] - gene[0] + 1}, start: {gene[0]}, end: {gene[1]}')
            for i, gene in enumerate(sequenceRegion)
        )
    # get sequence for intergenic region without concatenation
    else:
        sequence = (SeqRecord(
            Seq(seq[gene[0] - 1:gene[1]]),
            id="intergenic " + str(i + 1),
            description=f'length= {gene[1] - gene[0] + 1}, start: {gene[0]}, end: {gene[1]}')
            for i, gene in enumerate(sequenceRegion)
        )

    #inlist = list(sequence)
    #print(inlist[1])

    return sequence


# return a concatenated gene/ intergenic list ,i.e ['A','T','G','A'.........]
def merge(sequence):
    merged = []
    for i in range(len(sequence)):
        merged += (list(sequence[i].seq))
    return merged


# return the nucleotide frequency table as a dictionary, i.e. {'A':0.12, 'T':0.45,...}
def nuc_frenq(list_sequence, interLength):
    table = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for nuc in table:
        for i in list_sequence:
            if i == nuc:
                table[nuc] += 1
            else:
                table[nuc] += 0

        table[nuc] = table[nuc] / interLength
        table[nuc] = round(table[nuc], 4)
    return table


# return pandas dataframe of codon frequency of whole genic region (start+ middle+ stop)
# all rows printed (64 rows)
def total_cod(geneList, genLength):
    onestr = ''.join(geneList)
    amount = genLength / 3
    codons = (onestr[n:n + 3] for n in xrange(0, len(onestr), 3))  # creates generator
    dict_codons = {}
    for i in codons:
        if dict_codons.get(i) is not None:
            dict_codons[i] += 1
        else:
            dict_codons[i] = 1
    for i in dict_codons:
        dict_codons[i] /= amount
        dict_codons[i] = round(dict_codons[i], 3)
    df = pd.DataFrame(dict_codons.items(), columns=["codon", "frequency"])
    return df


# return pandas dataframe of codon frequency of start codon frequency
# default # of rows printed
def start_cod(listSequence, merged):
    seq = []
    for i in range(len(listSequence)):
        together = ''.join(str(listSequence[i].seq))
        seq += [together]
    startCodon = {}
    count = 0
    for i in seq:
        codons = ''.join((i[n:n + 3] for n in xrange(0, 3, 3)))
        count += 1
        if startCodon.get(codons) is None:
            startCodon.update({codons: 1})
        else:
            startCodon[codons] += 1
    for i in startCodon:
        startCodon[i] /= count
        startCodon[i] = round(startCodon[i], 3)
    onestr = ''.join(merged)
    codons = (onestr[n:n + 3] for n in xrange(0, len(onestr), 3))  # creates generator
    dict_codons = {}

    for codon in codons:
        if startCodon.get(codon) is not None:
            dict_codons[codon] = startCodon[codon]
        else:
            dict_codons[codon] = 0
    df2 = pd.DataFrame(dict_codons.items(), columns=["codon", "frequency"])
    return startCodon, df2


# return pandas dataframe of codon frequency of middle codon frequency
# default # of rows printed
def middle_cod(geneList, genLength):
    onestr = ''.join(geneList)
    amount = (genLength - 6) / 3
    codons = (onestr[n:n + 3] for n in xrange(3, len(onestr) - 3, 3))  # creates generator
    dict_codons = {}
    for i in codons:
        if dict_codons.get(i) is not None:
            dict_codons[i] += 1
        else:
            dict_codons[i] = 1
    # sum = 0
    for i in dict_codons:
        dict_codons[i] /= amount
        dict_codons['TAA'] = dict_codons['TGA'] = dict_codons['TAG'] = 0
        dict_codons[i] = round(dict_codons[i], 3)
        # sum += dict_codons[i]
    # print("total sum of frequency=: " + str(sumfre))
    pd.reset_option('display.max_rows')
    df = pd.DataFrame(dict_codons.items(), columns=["codon", "frequency"])
    return df


# return pandas dataframe of codon frequency of stop codon frequency
# default # of rows printed
def stop_cod(listSequence, merged):
    seq = []
    for i in range(len(listSequence)):
        together = ''.join(str(listSequence[i].seq))
        seq += [together]
    stopCodon = {}
    count = 0
    for i in seq:
        codons = ''.join((i[n:n + 3] for n in xrange(len(i) - 3, len(i), 3)))
        count += 1
        if stopCodon.get(codons) is None:
            stopCodon.update({codons: 1})
        else:
            stopCodon[codons] += 1
    for i in stopCodon:
        stopCodon[i] /= count
        stopCodon[i] = round(stopCodon[i], 3)
    onestr = ''.join(merged)
    codons = (onestr[n:n + 3] for n in xrange(0, len(onestr), 3))  # creates generator
    dict_codons = {}
    for codon in codons:
        if stopCodon.get(codon) is not None:
            dict_codons[codon] = stopCodon[codon]
        else:
            dict_codons[codon] = 0
    pd.reset_option('display.max_rows')
    df2 = pd.DataFrame(dict_codons.items(), columns=["codon", "frequency"])
    return stopCodon, df2


def file_genic(genic_sequence):
    SeqIO.write(genic_sequence, 'genic.fa', 'fasta')


def file_intergenic(intergenic_sequence):
    SeqIO.write(intergenic_sequence, 'intergenic.fa', 'fasta')


def main():
    a = Q1a()
    # 1
    print("1. Total genic length: " + str(a.totGenicLen))
    print("   Number of observed gene: " + str(len(a.genicRegion)))
    print("   Average genic length: " + str(a.avgenicLen) + "\n")

    # 2
    print("2. Total intergenic length: " + str(a.totInterLen))
    print("   Number of intergenic region: " + str(len(a.intergenicRegion)))
    print("   Average intergenic length= " + str(a.avgInterLen) + "\n")

    # 3
    print("3. The nucleotide frequency table is: " + str(a.nuc_frequency) + "\n")

    # 4
    pd.set_option('display.max_rows', None)
    print("\n" + "4. Genic codon frequency table as a whole: " + "\n" + str(a.totalCodon))

    # 5
    pd.reset_option('display.max_rows')
    print("\n" + "5. start codons are: " + str(a.startCodon))
    print("start Codon frequency table: " + "\n" + str(a.startTable))

    # 6
    pd.reset_option('display.max_rows')
    print("\n" + "6. Middle codon frequency table: " + "\n" + str(a.middleCodon))

    # 7
    pd.reset_option('display.max_rows')
    print("\n" + "7. stop codons are: " + str(a.stopCodon))
    print("stop Codon frequency table: " + "\n" + str(a.stopTable))


if __name__ == '__main__':
    main()
