import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def aln(str1, str2, read_name):
    align = pairwise2.align.localms(str1, str2, 2, -1, -1, -0.5)
    cnt = 0
    '''
    with open("Alignments.aln", 'a+') as file:
        file.write(read_name)
        file.write("\n")
        for a in align:
            file.write(format_alignment(*a))
    '''
    for letter in align[0][0]:
        if letter != '-':
            return cnt, int(align[0][2])
        cnt += 1

def alignment(str1, str2, read_name):
    if str1 > str2:
        return aln(str2, str1, read_name)
    return aln(str1, str2, read_name)


def get_score(str1, str2, read_name):
    '''
    aln = pairwise2.align.localms(str1, str2, 1, -1, -2, -0.5)
    
    with open("chunks.aln", 'a+') as outp:
        outp.write(read_name)
        outp.write("\n")
        for a in aln:
            outp.write(format_alignment(*a))
    '''
    return  int(pairwise2.align.localms(str1, str2, 1, -1, -2, -0.5, score_only=True)) #aln[0][2] #int(pairwise2.align.localms(str1, str2, 1, -1, -2, -0.5, score_only=True))



def main(str1, str2):
    aln = alignment(str1, str2)
    print(aln)


from Bio.pairwise2 import format_alignment


'''
def levenshtein(ref, seq):
    import numpy as np

    if len(ref) < len(seq):
        return levenshtein(seq, ref)
    if len(seq) == 0:
        return len(ref)

    ref = np.array(tuple(ref))
    seq = np.array(tuple(seq))

    previous_row = np.arange(seq.size + 1)
    for s in ref:
        current_row = previous_row + 1
#        print('s = %s\nPREVIOUS ROW = %s\nNEW ROW = %s' %(s, previous_row, current_row))
        current_row[1:] = np.minimum(current_row[1:], np.add(previous_row[:-1], seq != s))
        current_row[1:] = np.minimum(current_row[1:], current_row[0:-1] + 1)
        previous_row = current_row

    return previous_row[-1] #- (len(ref) - len(seq))


def levenshtein2(str1, str2):

    if len(str1) < len(str2):
        return levenshtein2(str2, str1)
    if len(str1) == 0:
        return len(str1)

    n = len(str1)
    p = len(str2)

    # создаем матрицу с оптимальными комбинациями пар
    mat = [[0 for i in range(n + 1)] for j in range(p + 1)]


    for i in range(len(str1) + 1):
        for j in range(len(str2) + 1):
            if i > 0 and j > 0:
                if str1[i - 1] == str2[j - 1]:
                    mat[i][j] = max([mat[i - 1][j], mat[i][j - 1], mat[i - 1][j - 1] + 1])
                else:
                    mat[i][j] = max([mat[i - 1][j], mat[i][j - 1], mat[i - 1][j - 1]])

    return mat[-1][-1]


def levenshtein3(s1, s2):
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1  # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1  # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1] - (len(s1) - len(s2))


def edDistDp(x,	y):
    """	Calculate edit distance between	sequences x and	y using matrix dynamic programming. Return	distance.	"""
    D = np.zeros((len(x)+1,	len(y)+1),	dtype=int)
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)
    return D[len(x), len(y)] - (len(y) - len(x))


def iterative_levenshtein(s, t):
    """ 
        iterative_levenshtein(s, t) -> ldist
        ldist is the Levenshtein distance between the strings 
        s and t.
        For all i and j, dist[i,j] will contain the Levenshtein 
        distance between the first i characters of s and the 
        first j characters of t
    """
    rows = len(s) + 1
    cols = len(t) + 1
    dist = [[0 for x in range(cols)] for x in range(rows)]
    # source prefixes can be transformed into empty strings
    # by deletions:
    for i in range(1, rows):
        dist[i][0] = i
    # target prefixes can be created from an empty source string
    # by inserting the characters
    for i in range(1, cols):
        dist[0][i] = i

    for col in range(1, cols):
        for row in range(1, rows):
            if s[row - 1] == t[col - 1]:
                cost = 0
            else:
                cost = 1
            dist[row][col] = min(dist[row - 1][col] + 1,  # deletion
                                 dist[row][col - 1] + 1,  # insertion
                                 dist[row - 1][col - 1] + cost)  # substitution
    #for r in range(rows):
        #print(dist[r])

    return dist[row][col] #- len(t) + len(s)
'''



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