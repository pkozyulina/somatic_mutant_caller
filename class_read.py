import collections
import re
from Bio import SeqIO
from levenshtein import alignment
from levenshtein import get_score

class Read:

    def __init__(self, seq, name, quality, strand):
        self.seq = seq
        self.name = name  # name of the read
        self.norm = []  # self.norm[0] - name of a ref contig that read has aligned to
                        # self.norm[1] - relative position on the reference
                        # self.norm[2] - alignment score
        self.mut = collections.defaultdict(tuple)  # key = mutation, value = (position, alignment score)
        self.quality = quality
        self.coordinates = None  # global coordinates of the read
        self.possible_mutations = []  # list of tuples: (mutations that are covered by read, positions of the mutations)
        self.strand = strand  # 0 - forward strand, 1 - reverse strand

    def quant_dist_norm(self, ref_candidates, ref_dict):
        min_dist = -1000000000
        ref_min = ''
        pos_min = None

        if not ref_candidates:
            ref_candidates = ref_dict

        for ref in ref_candidates:
            pos, dist = alignment(self.seq, ref_dict[ref][0], self.name)
            if dist > min_dist:
                min_dist = dist
                ref_min = ref
                pos_min = pos

        self.norm = [ref_min, pos_min, min_dist]

        return ref_min


    def quant_dist_mut(self, ref_dict):
        for mutation, _ in self.possible_mutations:
            self.mut[mutation] = alignment(self.seq, ref_dict[mutation][0], self.name)


    def is_mut(self):
        for mutation, _ in self.possible_mutations:
            if self.mut[mutation][1] > self.norm[2]:
                return mutation
        return None


    def calc_delta(self, mut_pos_global):
        delta_start = mut_pos_global - self.coordinates[0]

        if delta_start >= 5 and delta_start <= len(self.seq) - 5:
            delta_qual = self.quality[delta_start - 5:delta_start + 5]
        elif delta_start < 5:
            delta_qual = self.quality[0:delta_start + 5]
        elif delta_start > len(self.seq) - 5:
            delta_qual = self.quality[delta_start - 5:]

        return sum(delta_qual)/len(delta_qual)


    # define global coordinates for read
    def set_coordinates(self, reference):
        self.coordinates = reference.position_norm(self.norm[0])
        self.coordinates[0] += self.norm[1]
        self.coordinates[1] = self.coordinates[0] + len(self.seq)


    # check if read overlaps any of mutations
    def is_correct_pos(self, reference):

        self.set_coordinates(reference)

        for mutation, pos in reference.position_mut(self.norm[0]):

            if pos >= self.coordinates[0] and pos <= self.coordinates[1]:
                self.possible_mutations.append((mutation, pos))

        if self.possible_mutations:
            return True

        return False


    def __str__(self):
        return str(self.seq)

    def __len__(self):
        return len(self.seq)




# Creating a reference dictionary
class Reference:

    def __init__(self):
        self.norm = collections.defaultdict(tuple) # key = reference name; value = (reference sequence,
                                                   # reference global position), where reference global position is a tuple
        self.mutant = collections.defaultdict(dict) # key = reference name; value = {}, where key = mutation, value =
                                                    # (reference sequence, reference global position, allele lengths ),
                                                    # where allele lengths is a tuple


    # getting global coordinates of each reference
    def get_position(self, line):
        pattern1 = r'([^|]*)'
        ref_name = (re.match(pattern1, line)).group(1)
        status = False

        if "|mut|" in line:  # here we check if its mutant or not
            pattern2 = r'(\|mut\|[^|]*\|[^|]*)\|(\d+)\|'
            tmp = re.search(pattern2, line)
            ref_position = int(tmp.group(2))
            status = True

        if "|wt|" in line:
            pattern2 = r'\|(\d+)\|(\d+)$'
            tmp = re.search(pattern2, line)
            ref_position = (int(tmp.group(1)), int(tmp.group(2)))

        return status, ref_name, ref_position


    # adding length of normal and mutant alleles
    def allele_length(self, line):

        pattern = r'\|([a-zA-Z-]+)\|([a-zA-Z-]+)$'
        tmp = re.search(pattern, line)

        norm_len = len(tmp.group(1))
        mut_len = len(tmp.group(2))

        if tmp.group(1) == "-":
            mut_len += 1
        if tmp.group(2) == "-":
            norm_len += 1

        return (norm_len, mut_len)


    # parsing reference file
    def parse_reference(self, ref_file):
        for rec in SeqIO.parse(ref_file, "fasta"):
            line = rec.id
            status, ref_name, ref_position = self.get_position(line)

            if status:
                alleles_len = self.allele_length(line)
                self.mutant[ref_name][rec.id] = (rec.seq, ref_position, alleles_len)

            else:
                self.norm[ref_name] = (rec.seq, ref_position)


    # looking for a suitable candidates to map to
    def try_ref(self, chunk):

        dist_max = -999999
        dist_all = []
        all_max = []

        for ref in self.norm:
            dist = get_score(chunk, self.norm[ref][0], chunk)
            dist_all.append((ref, dist))
            if dist > dist_max:
                dist_max = dist

        for ref, dist in dist_all:
            if dist == dist_max:
                all_max.append(ref)

        return all_max, dist_max


    # global coordinates of normal "gene"
    def position_norm(self, ref_name):
        return [self.norm[ref_name][1][0], self.norm[ref_name][1][1]]


    # global coordinate of mutation
    def position_mut(self, ref_name):
        for key, value in self.mutant[ref_name].items():
            yield key, value[1]

    # calculating mutant allele length for probability quantification
    def mutation_length(self, ref_name, mutant):
        n, m = self.mutant[ref_name][mutant][2]
        return n, m