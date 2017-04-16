import argparse
import collections
import re
from levenshtein import alignment
from class_read import Read
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
            yield record.seq, record.letter_annotations["phred_quality"]


# Creating a reference dictionary
def parse_reference(ref_file):
    ref_norm_dict = collections.defaultdict(dict)
    ref_mut_dict = collections.defaultdict(dict)
    for rec in SeqIO.parse(ref_file, "fasta"):
        line = rec.id
        if "|" in line:  # here we check if its mutant or not
            pattern = r'([^|]*)'
            ref_name = (re.match(pattern, line)).group(1)
            ref_mut_dict[ref_name][rec.id] = rec.seq
        else:
            ref_name = rec.id
            ref_norm_dict[ref_name] = rec.seq
    return ref_norm_dict, ref_mut_dict


# quantifing distances for norm and mut
def dist_quant(seq, quality, ref_candidates, ref_norm_dict, ref_mut_dict):
    read = Read(seq, quality)
    ref = read.quant_dist_norm(ref_candidates, ref_norm_dict)
    if not ref:
        return None, None
    read.quant_dist_mut(ref_mut_dict[ref])
    return ref, read


# looking for a suitable candidates to map to
def try_ref(seq, ref_dict):
    ref_min = []
    for ref in ref_dict:
        dist = alignment(ref_dict[ref], seq)
        if dist >= -4:  # this number is still needed to be adjusted
            ref_min.append(ref)
    return ref_min


# here we will have a dict with gene names and how many reads aligns to norm and mutant form
def build_alignment(input_file, ref_file, qual):

    gene_list = collections.defaultdict(list)
    ref_norm_dict, ref_mut_dict = parse_reference(ref_file)

    for seq, quality in tqdm(parse_fastq(input_file, qual)):

        # roughly checking where the read maps
        chunk = seq[round(len(seq) / 2 - 7):round(len(seq) / 2 + 7)]
        ref_candidates = try_ref(chunk, ref_norm_dict)
        ref_candidates.extend(try_ref(twin(chunk), ref_norm_dict))

        # getting an actual reference mapped
        ref, read = dist_quant(seq, quality, ref_candidates, ref_norm_dict, ref_mut_dict)
        ref_r, read_r = dist_quant(twin(seq), quality, ref_candidates, ref_norm_dict, ref_mut_dict)

        if ref:
            gene_list[ref].append(read)
        if ref_r:
            gene_list[ref_r].append(read_r)

    return gene_list  # returns a dict with all contigs as keys and a list of aligned reads to each contig




'''Object of class Read has info about levenshtein distance to normal and all similar mutant references. 
This function gets all Reads that align to the same location ("gene") and count number of reads aligned to norm and to mut references. 
We feed this function with read list for each contig from build alignment function (gene_list[ref])'''

# counting number of mutant and normal reads mapped to the reference
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
    mu = 10 ** ((mu / i) / -10)  # lambda of poisson distribution
    return mu, x


def poisson_prob(mu, x):
    return poisson.cdf(k=x, mu=mu)  # calculating cumulative distribution function (probability?)


def main():

    # parsing command line arguments
    parser = argparse.ArgumentParser(description='Clinically relevant somatic mutations caller tool')
    parser.add_argument('-i', '--input', help='Input fastq file', metavar='File', type=argparse.FileType(),
                        required=True)
    parser.add_argument('-r', '--reference', help='Amplicon reference in fasta-format for both normal and mutant forms',
                        metavar='File',
                        type=argparse.FileType(), required=True)
    parser.add_argument('-q', '--quality', help='Quality filtering (default: 30)', metavar='Int', type=int,
                        default=30)
    parser.add_argument('-o', '--output', help='Output vcf', metavar='File', type=argparse.FileType('w'),
                        required=True)
    args = parser.parse_args()


    # reading file and bulding a dict with all the distances
    align = build_alignment(args.input, args.reference, args.quality)


    # running the stats
    mutant_vs_norm_dict = collections.defaultdict(list)

    for ref, reads in align.items():
        mu, x = mut_vs_norm(reads)

        for mutant in x:
            if mutant != ref:
                prob = poisson_prob(mu, x[mutant])
                mutant_vs_norm_dict[mutant] = [x[mutant], sum(x.values()),
                                               prob]  # mutation, mutation counts (number of reads), sum counts for the reference


    # writing a table with Levenshtein distances to norm and mut references per read
    with open(args.output.name, 'w') as outpa:
        i = 0
        outpa.write('Read\tAligned to\tNorm dist\tMut\tMut dist\tis mutant?\n')
        for key, values in align.items():
            for value in values:
                i += 1
                for_print = ''
                mutant = value.is_mut()

                if not mutant:
                    mutant = "Norm"

                for mut, dist in value.mut.items():
                    for_print += '%s\t%i\t' % (mut, dist)
                outpa.write('%i\t%s\t%i\t%s\t%s\n' % (i, value.norm[0], value.norm[1], for_print, mutant))


    # writing a table with Poisson probability for mutations
    with open("4005_probability_table.tsv", 'w') as outpb:
        i = 0
        outpb.write('N\tMutation\tMut counts\tSum counts\tProbability\n')
        for key, values in mutant_vs_norm_dict.items():
            i += 1
            outpb.write('%i\t%s\t%i\t%i\t%s\n' % (i, key, values[0], values[1], values[2]))


if __name__ == '__main__':
    main()
