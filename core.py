import argparse
import collections
import numpy as np
from levenshtein import levenshtein
from scipy.stats import poisson
from Bio import SeqIO
from Bio import Seq
import re
from tqdm import tqdm




class Read:

    def __init__(self, seq, quality):
        self.seq = seq
        self.norm = []  # self.norm[0] - name of a ref contig that read has aligned to and self.norm[1] - Levenshtein distance to it
        self.mut = collections.defaultdict(int)  # all distances to mutant contigs
        self.quality = quality
        self.avg_quality = sum(quality)/len(quality)


    def quant_dist_norm(self, ref_dict):
        min_dist = 1000
        ref_min = ''

        for ref in ref_dict:
            dist = levenshtein(ref_dict[ref], self.seq)
            if dist < min_dist:
                min_dist = dist
                ref_min = ref

        self.norm = [ref_min, min_dist]

        return ref_min


    def quant_dist_mut(self, ref_dict):
        for ref in ref_dict:
            self.mut[ref] = levenshtein(ref_dict[ref], self.seq)


    def is_mut(self):
        for key, value in self.mut:
            if value < self.norm[2]:
                return key
        return None

    def __str__(self):
        return str(self.seq)

    def __len__(self):
        return len(self.seq)





# creating reverse complement read
def twin(read):
    return Seq.reverse_complement(read)


# Parsing fastq
def parse_fastq(input_file, quality):
    for record in SeqIO.parse(input_file, "fastq"):
        if min(record.letter_annotations["phred_quality"]) >= quality:
            yield record.seq, record.letter_annotations["phred_quality"]


# Creating a reference dictionary
def parse_reference(ref_file):
    record_dict = collections.defaultdict(dict)
    for rec in SeqIO.parse(ref_file, "fasta"):
        line = rec.id
        if "|" in line:
            pattern = r'([^|]*)'
            ref_name = (re.match(pattern, line)).group(1)
            record_dict[ref_name][rec.id] = rec.seq
        else:
            ref_name = rec.id
            record_dict[ref_name] = rec.seq
    return record_dict


# quantifing distances for norm and mut
def dist_quant(seq, quality, ref_norm_dict, ref_mut_dict):
    read = Read(seq, quality)
    ref = read.quant_dist_norm(ref_norm_dict)
    read.quant_dist_mut(ref_mut_dict[ref])
    return ref, read


# here we will have a dict with gene names and how many reads aligns to norm and mutant form
def build_alignment(input_file, ref_file_norm, ref_file_mut, quality):

    gene_list = collections.defaultdict(list)
    ref_norm_dict = parse_reference(ref_file_norm)
    ref_mut_dict = parse_reference(ref_file_mut)

    for seq, quality in tqdm(parse_fastq(input_file, quality)):
        ref, read = dist_quant(seq, quality, ref_norm_dict, ref_mut_dict)
        ref_r, read_r = dist_quant(twin(seq), quality, ref_norm_dict, ref_mut_dict)
        gene_list[ref].append(read)
        gene_list[ref_r].append(read_r)

    return gene_list  # returns a dict with all contigs as keys and a list of aligned reads to each contig




'''Object of class Read has info about levenshtein distance to normal and all similar mutant references. 
This function gets all Reads that align to the same location ("gene") and count number of reads aligned to norm and to mut references. 
We feed this fucntion with read list for each contig from build alignment function (gene_list[ref])'''

def mut_vs_norm(reads):
    mu, i = 0, 0
    x = collections.defaultdict(int)  # x[norm] - normal read counts and x[mut] - mutant read counts
    for read in reads:
        i += 1
        mu += read.avg_quality  # calculating sum of read average quality
        mutant = read.is_mut()
        if mutant:
            x[mutant] += 1  # mutant read counts
        else:
            x[read.norm[0]] += 1  # normal read counts
    mu = (10 ** (mu / i)) / (-10)  # lambda of poisson distribution
    return mu, x


def poisson_prob(mu, x, size):
    return poisson.cdf(k=x, mu=mu, size=size)  # calculating cumulative distribution function





def main():

    # parsing command line arguments
    parser = argparse.ArgumentParser(description='Clinically relevant somatic mutations caller tool')
    parser.add_argument('-i', '--input', help='Input fastq file', metavar='File', type=argparse.FileType(),
                                    required=True)
    parser.add_argument('-n', '--nref', help='Normal amplicon reference in fasta-format', metavar='File',
                                    type=argparse.FileType(), required=True)
    parser.add_argument('-m', '--mref', help='Mutant amplicon reference in fasta-format', metavar='File',
                                    type=argparse.FileType(), required=True)
    parser.add_argument('-q', '--quality', help='Quality filtering (default: 30)', metavar='Int', type=int,
                                    default=30)
    parser.add_argument('-o', '--output', help='Output vcf', metavar='File', type=argparse.FileType('w'),
                                    required=True)
    args = parser.parse_args()


    # reading file and bulding a dict with all the distances
    align = build_alignment(args.input, args.nref, args.mref, args.quality)
    probability_dict = collections.defaultdict(int)
    for ref, reads in align.items():
        mu, x = mut_vs_norm(reads)
        for mutant in x:
            if mutant != ref:
                prob = poisson_prob(mu, x[mutant], sum(x.values()))
                probability_dict[mutant] = prob

    with open(args.output.name, 'w') as outp:
        i = 0
        outp.write('Read\tAligned to\tNorm dist\tMut\tMut dist\n')
        for key, values in align.items():
            for value in values:
                i += 1
                for_print = ''
                for mut, dist in value.mut.items():

                    for_print += '%s\t%i\t' % (mut, dist)
                outp.write('%i\t%s\t%i\t%s\n' % (i, value.norm[0], value.norm[1], for_print))


if __name__ == '__main__':
    main()