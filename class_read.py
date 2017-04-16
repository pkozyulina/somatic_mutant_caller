import collections
from levenshtein import alignment

class Read:

    def __init__(self, seq, quality):
        self.seq = seq
        self.norm = []  # self.norm[0] - name of a ref contig that read has aligned to and self.norm[1] - Levenshtein distance to it
        self.mut = collections.defaultdict(int)  # all distances to mutant contigs
        self.quality = quality
        self.avg_quality = sum(quality)/len(quality)


    def quant_dist_norm(self, ref_candidates, ref_dict):
        min_dist = -1000000000000
        ref_min = ''

        if not ref_candidates:
            ref_candidates = ref_dict

        for ref in ref_candidates:
            dist = alignment(self.seq, ref_dict[ref])
            if dist > min_dist:
                min_dist = dist
                ref_min = ref

        self.norm = [ref_min, min_dist]

        if abs(min_dist) > 30:
            return None

        return ref_min


    def quant_dist_mut(self, ref_dict):
        for ref in ref_dict:
            self.mut[ref] = alignment(self.seq, ref_dict[ref])


    def is_mut(self):
        for key, value in self.mut.items():
            if value > self.norm[1]:
                return key
        return None


    def __str__(self):
        return str(self.seq)


    def __len__(self):
        return len(self.seq)
