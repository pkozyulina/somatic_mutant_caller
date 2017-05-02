import collections
import re
from Bio import SeqIO
from levenshtein import alignment
from levenshtein import get_score

class Read:

    def __init__(self, seq, quality):
        self.seq = seq
        self.norm = []  # self.norm[0] - name of a ref contig that read has aligned to
                        # self.norm[1] - relative position on the reference
                        # self.norm[2] - alignment score
        self.mut = collections.defaultdict(tuple)  # key = mutation, value = (position, alignment score)
        self.quality = quality
        self.coordinates = None  # global coordinates of the read
        self.possible_mutations = []  # list of tuples: (mutations that are covered by read, positions of the mutations)


    def quant_dist_norm(self, ref_candidates, ref_dict):
        min_dist = -1000000000
        ref_min = ''
        pos_min = None

        if not ref_candidates:
            ref_candidates = ref_dict

        for ref in ref_candidates:
            pos, dist = alignment(self.seq, ref_dict[ref][0])
            if abs(dist) < abs(min_dist):
                min_dist = dist
                ref_min = ref
                pos_min = pos

        self.norm = [ref_min, pos_min, min_dist]

        return ref_min


    def quant_dist_mut(self, ref_dict):
        for mutation, _ in self.possible_mutations:
            self.mut[mutation] = alignment(self.seq, ref_dict[mutation][0])


    def is_mut(self):
        for mutation, _ in self.possible_mutations:
            if abs(self.mut[mutation][1]) < abs(self.norm[2]):
                return mutation
        return None


    def calc_delta(self, mut_pos_global):
        delta_start = mut_pos_global - self.coordinates[0]
        delta_qual = self.quality[delta_start - 5:delta_start + 5]
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
            if pos > self.coordinates[0] and pos < self.coordinates[1]:
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
                                                    # (reference sequence, reference global position), where reference global position is a tuple

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


    def parse_reference(self, ref_file):
        for rec in SeqIO.parse(ref_file, "fasta"):
            line = rec.id
            status, ref_name, ref_position = self.get_position(line)
            if status:
                self.mutant[ref_name][rec.id] = (rec.seq, ref_position)
            else:
                self.norm[ref_name] = (rec.seq, ref_position)


    # looking for a suitable candidates to map to
    def try_ref(self, chunk):

        dist_min = 999999
        dist_all = []
        all_min = []

        for ref in self.norm:
            dist = get_score(chunk, self.norm[ref][0])
            dist_all.append((ref, dist))
            if abs(dist) < abs(dist_min):
                dist_min = dist

        for ref, dist in dist_all:
            if dist == dist_min:
                all_min.append(ref)

        return all_min, dist_min


    def position_norm(self, ref_name):
        return [self.norm[ref_name][1][0], self.norm[ref_name][1][1]]


    def position_mut(self, ref_name):
        for key, value in self.mutant[ref_name].items():
            yield key, value[1]