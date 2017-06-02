#!/usr/bin/env python3

import argparse
import collections
import re
from levenshtein import alignment
from levenshtein import get_score
from class_read import Read
from class_read import Reference
from scipy.stats import poisson
from Bio import SeqIO
from Bio import Seq
from tqdm import tqdm


# creating reverse complement read
def twin(read):
    return Seq.reverse_complement(read)


# Parsing fastq
def parse_fastq(input_file, quality):
    for record in SeqIO.parse(input_file, "fastq"):
        if min(record.letter_annotations["phred_quality"]) >= quality:
            yield str(record.seq), record.letter_annotations["phred_quality"], record.id


# quantifing distances for norm and mut
def dist_quant(seq, quality, name, ref_candidates, reference, strand):
    read = Read(seq, name, quality, strand)
    ref = read.quant_dist_norm(ref_candidates, reference.norm)
    if read.is_correct_pos(reference):
        read.quant_dist_mut(reference.mutant[ref])
        return ref, read
    return None, None

def check_chunk(seq, reference):

    idx = int(len(seq)/2)
    # roughly checking where the read maps
    chunk = seq[idx - 15:idx + 16]
    ref_candidate, dist_fw = reference.try_ref(chunk)
    ref_candidate_r, dist_rv = reference.try_ref(twin(chunk))

    score_threshold = int(len(chunk)*0.6)

    if dist_rv < score_threshold and dist_fw < score_threshold:
        return None, None

    # choosing a strand of read (forward or revers)
    '''
    if dist_fw == dist_rv:
        chunk = seq
        ref_candidate, dist_fw = reference.try_ref(chunk)
        ref_candidate_r, dist_rv = reference.try_ref(twin(chunk))
    '''
    if dist_fw < dist_rv:
        return ref_candidate_r, 1

    return ref_candidate, 0  # here we roughly decide which reference read aligns to



# here we will have a dict with gene names and how many reads aligns to norm and mutant form
def build_alignment(input_files, reference, qual):

    gene_list = collections.defaultdict(list)
    for file in input_files:
        for seq, quality, name in tqdm(parse_fastq(file, qual)):
            ref_candidate, strand = check_chunk(seq, reference)
            if ref_candidate:
                if strand:
                    ref, read = dist_quant(twin(seq), quality, name, ref_candidate, reference, 1)  # getting an actual reference mapped
                else:
                    ref, read = dist_quant(seq, quality,  name, ref_candidate, reference, 0)

                if ref:
                    gene_list[ref].append(read)

    return gene_list  # returns a dict with all contigs as keys and a list of aligned reads to each contig


'''Object of class Read has info about levenshtein distance to normal and all similar mutant references. 
This function gets all Reads that align to the same location ("gene") and count number of reads aligned to norm and to mut references. 
We feed this function with read list for each contig from build alignment function (gene_list[ref])'''

def collect_mutants(reads, ref, reference):

    mu = collections.defaultdict(int)
    x = collections.defaultdict(int)  # x[norm] - normal read counts and x[mut] - mutant read counts
    cnt = collections.defaultdict(int)

    # calculating a number of mutant reads vs normal ones
    for read in reads:
        if read.is_correct_pos(reference):

            for mutation, position in read.possible_mutations:

                mu[mutation] += read.calc_delta(position)  # calculating sum of delta average qualities for every possible mutation position
                cnt[mutation] += 1

            mutant = read.is_mut()
            if mutant:
                x[mutant] += 1  # mutant read counts
            else:
                x[read.norm[0]] += 1  # normal read counts

    return mu, x, cnt


def mut_vs_norm(reads, ref, reference):

    mu, x, cnt = collect_mutants(reads, ref, reference)

    # quantifing lambda poisson distr for each
    to_delete = []
    for mutation in mu:
        if mutation in x:
            mu[mutation] = 10 ** ((mu[mutation] / cnt[mutation]) / -10)  # lambda of poisson distribution
        else:
            to_delete.append(mutation)

    # deleting mutations that are not present
    for mutation in to_delete:
        del(mu[mutation])

    return mu, x


def poisson_prob(mu, x, n, m):  ## n = length of reference allele, m = length of mutant allele
    return ((1 - poisson.cdf(k=x, mu=mu))**n)/(4**m - 1)  # calculating cumulative distribution function (probability)


# running the stats
# actually creates all objects and do all stats
def stats(input_files, ref_file, input_quality):

    # create reference
    reference = Reference()
    reference.parse_reference(ref_file)

    # reading file and bulding a dict with all the distances
    align = build_alignment(input_files, reference, input_quality)
    mutant_vs_norm_dict = collections.defaultdict(list)

    for ref, reads in align.items():
        mu, x = mut_vs_norm(reads, ref, reference)

        for mutant in x:
            if mutant != ref:
                n, m = reference.mutation_length(ref, mutant)
                prob = poisson_prob(mu[mutant], x[mutant], n, m)
                mutant_vs_norm_dict[mutant] = [x[mutant], sum(x.values()), prob]
                # mutation, mutation counts (number of reads), sum counts for the reference

    return align, mutant_vs_norm_dict


def main():

    # parsing command line arguments
    parser = argparse.ArgumentParser(description='Clinically relevant somatic mutations caller tool')
    parser.add_argument('-i', '--input', help='Input fastq file', metavar='File', type=argparse.FileType(), nargs='+',
                        required=True)
    parser.add_argument('-r', '--reference', help='Amplicon reference in fasta-format for both normal and mutant forms',
                        metavar='File',
                        type=argparse.FileType(), required=True)
    parser.add_argument('-q', '--quality', help='Quality filtering (default: 30)', metavar='Int', type=int,
                        default=30)
    parser.add_argument('-o', '--output', help='Two output names: for table with reads and table with probablities', metavar='File', type=argparse.FileType('w'), nargs=2,
                        required=True)
    args = parser.parse_args()



    # here we invoke all the functions needed
    align, mutant_vs_norm_dict = stats(args.input, args.reference, args.quality)


    # writing a table with Alignment scores to norm and mut references per read
    with open(args.output[0].name, 'w') as outpa:
        outpa.write('Read\tStrand\tAligned to\tNorm dist\tMut\tMut dist\tis mutant?\n')
        for key, values in align.items():
            for value in values:
                mutant = value.is_mut()

                if not mutant:
                    mutant = "Norm"

                outpa.write('%s\t%s\n' % (value.print_all(), mutant))


    # writing a table with Poisson probability for mutations
    with open(args.output[1].name, 'w') as outpb:
        i = 0
        outpb.write('N\tMutation\tMut counts\tSum counts\tProbability\n')
        for key, values in mutant_vs_norm_dict.items():
            i += 1
            outpb.write('%i\t%s\t%i\t%i\t%s\n' % (i, key, values[0], values[1], values[2]))


if __name__ == '__main__':
    main()