
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

    return previous_row[-1] - (len(ref) - len(seq))

def levenshtein2(str1, str2):

    if len(str1) < len(str2):
        return levenshtein2(str2, str1)
    if len(str1) == 0:
        return len(str1)

    n = len(str1)
    p = len(str2)
    print(n, p)

    # создаем матрицу с оптимальными комбинациями пар
    mat = [[0 for i in range(n + 1)] for j in range(p + 1)]
    print(mat)

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

def alignment(str1, str2):
    if len(str1) < len(str2):
        return alignment(str2, str1)
    if len(str1) == 0:
        return len(str1)

    n = len(str1)
    p = len(str2)

    # создаем матрицу с оптимальными комбинациями пар
    mat = [[0 for i in range(n + 1)] for j in range(p + 1)]
    idx = (0, 0)
    maxi = 0
    for i in range(len(str1) + 1):
        for j in range(len(str2) + 1):
            if i > 0 and j > 0:
                if str1[i - 1] == str2[j - 1]:
                    mat[i][j] = max([mat[i - 1][j], mat[i][j - 1], mat[i - 1][j - 1] + 1])
                else:
                    mat[i][j] = max([mat[i - 1][j], mat[i][j - 1], mat[i - 1][j - 1]])
                if mat[i][j] > maxi:
                    maxi = mat[i][j]
                    idx = [i, j]

    # создаем "схему" будущего выравнивания, где 0 = совпадение, 1 = вставка, -1 = делеция
    path = []
    i, j = idx
    while i != 0 or j != 0:
        if i == 0:
            path.append(1)
            j -= 1
            continue
        elif j == 0:
            path.append(-1)
            i -= 1
            continue

        if str1[i-1] == str2[j-1]:
            path.append(0)
            i -= 1
            j -= 1
        elif mat[i-1][j] > mat[i][j-1]:
            path.append(-1)
            i -= 1
        else:
            path.append(1)
            j -= 1

    # разворачиваем путь спереди назад (чтобы удобнее было бежать)
    path.reverse()
    ans = {str1: '', str2: ''}
    k = 0
    l = 0
    # бежим по циклу и расшифровываем наш "код"
    for el in path:
        if el == 0:
            ans[str1] += str1[k]
            ans[str2] += str2[l]
            k += 1
            l += 1
        elif el == 1:
            ans[str1] += '-'
            ans[str2] += str2[l]
            l += 1
        elif el == -1:
            ans[str1] += str1[k]
            ans[str2] += '-'
            k += 1

    return mat[-1][-1], ans[str1], ans[str2] # возвращает количество соответствий, а также первую и вторую строку выравнивания


def main(str1, str2):
    print(levenshtein(str1, str2))

if __name__ == '__main__':
    import sys
    str1 = "TTTTCGTGGAAGTGGGTTA" #CCTGACAGTGTGCACGCCCCCAGCAGGTTCACAATATTCTCGTGGCTTCCCAGCTGGGTCATCATCTTGAGTTCTGACATGAGTGCCTCTCTTTCAGAGCTGTCTGCTTTTTCTGTCAAAGAAAGGAGCATT"
    #str2 = "CTTCCCCCGCCCCCGCAGTCTCTCTCTAACTCTCTCTCCTTCCCCCAACCCTCCC"
    file = {"ABL1_1": "GGTGTGTCCCCCAACTACGACAAGTGGGAGATGGAACGCACGGACATCACCATGAAGCACAAGCTGGGCGGGGGCCAGTACGGGGAGGTGTACGAGGGCGTGTGGAAGAAATACAGCCTGACGGTGGCCGTGAAGACCTTGAAGGTAGGCTGGGACTGCCGGGGGTGCCCAGGGTACGTGGGGCAAGGCGTCTGCTGGCATTAGGCGATGCATCTGCCTGGAAGTCTACCTCCTGCC",
            "KDR6": "TGCCATAGCATGCAGGAAGCACTAGCCAGTACCTTCCTCTTCTTCTACATCACTGAGGGACTTCTCCTCCACAAATCCAGAGCTGGCTGAGCTCTGGCTACTGGTGATGCTGTCCAAGCGCCGTTTCAGATCCACAGGGATTGCTCCAACGTAGTCTTTCCCTTGACGGAATCGTGCCCCTTTGGTCTATAAAAAAGCAAAGGAACAAACAAACTC",
            "TP53_1": "CTCCTTCCCAGCCTGGGCATCCTTGAGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGCAGCAACAGAGGAGGGGGAGAAGTAAGTATATACACAGTACCTGAGTTAAAAGATGGTTCAAGTTACAATTGTTTGACTTTATGACGGTACAAAAGCAACATGCATTTAGTAGAAACTGCACTTCAAGTACCTA"}
    for key, str2 in file.items():
        print(key)
        main(str1, str2)