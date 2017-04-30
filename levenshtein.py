import numpy as np
from Bio import pairwise2


def aln(str1, str2):
    align = pairwise2.align.localms(str1, str2, 1, -1, -2, -1)
    cnt = 0
    for letter in align[0][0]:
        if letter != '-':
            return cnt, int(align[0][2]) # - len(str1)
        cnt += 1

def alignment(str1, str2):
    if str1 > str2:
        return aln(str2, str1)
    return aln(str1, str2)


def get_score(str1, str2):
    return int(pairwise2.align.localms(str1, str2, 1, -1, -2, -1, score_only=True))# - len(str1))



def main(str1, str2):
    aln = alignment(str1, str2)
    print(aln)
    #print(aln[10])

from Bio.pairwise2 import format_alignment



if __name__ == '__main__':

    import sys
    str1 = "TTGAGTTCCAAGGCCTCTTCTTCTACATCACTGAGGGACTT"
    file = {"ABL1_1": "GGTGTGTCCCCCAACTACGACAAGTGGGAGATGGAACGCACGGACATCACCATGAAGCACAAGCTGGGCGGGGGCCAGTACGGGGAGGTGTACGAGGGCGTGTGGAAGAAATACAGCCTGACGGTGGCCGTGAAGACCTTGAAGGTAGGCTGGGACTGCCGGGGGTGCCCAGGGTACGTGGGGCAAGGCGTCTGCTGGCATTAGGCGATGCATCTG",
            "KDR6":   "TGCCATAGCATGCAGGAAGCACTAGCCAGTACCTTCCTCTTCTTCTACATCACTGAGGGACTTCTCCTCCACAAATCCAGAGCTGGCTGAGCTCTGGCTACTGGTGATGCTGTCCAAGCGCCGTTTCAGATCCACAGGGATTGCTCCAACGTAGTCTTTCCCTTGACGGAATCGTGCCCCTTTGGTCTATAAAAAAGCAAAGGAACAAACAAACTC",
            "TP53_1": "CTCCTTCCCAGCCTGGGCATCCTTGAGTTCCAAGGCCTCTTCTTCTACATCACTGAGGGACTTAGCGCTCACGCCCACGGATCTGCAGCAACAGAGGAGGGGGAGAAGTAAGTATATACACAGTACCTGAGTTAAAAGATGGTTCAAGTTACAATTGTTTGACTTTATGACGGTACAAAAGCAACATGCATTTAGTAGAAACTGCACTTCAAGTACC"}

    for key, str2 in file.items():
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n', key)
        #print(len(str1), len(str2))
        main(str1, str2)
