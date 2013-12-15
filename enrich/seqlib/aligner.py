from __future__ import print_function
import numpy as np


_simple_similarity = {
        'A' : {'A' : 1, 'C' : -1, 'G' : -1, 'T' : -1, 'N' : 0, 'X' : 0},
        'C' : {'A' : -1, 'C' : 1, 'G' : -1, 'T' : -1, 'N' : 0, 'X' : 0},
        'G' : {'A' : -1, 'C' : -1, 'G' : 1, 'T' : -1, 'N' : 0, 'X' : 0},
        'T' : {'A' : -1, 'C' : -1, 'G' : -1, 'T' : 1, 'N' : 0, 'X' : 0},
        'N' : {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0, 'N' : 0, 'X' : 0},
        'X' : {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0, 'N' : 0, 'X' : 0},
        'gap' : -1
}


class Aligner(object):
    _MAT = 1    # match
    _INS = 2    # insertion (with respect to wild type)
    _DEL = 3    # deletion (with respect to wild type)
    _END = 4    # end of traceback

    def __init__(self, similarity=_simple_similarity):
        similarity_keys = similarity.keys()
        if 'gap' in similarity_keys:
            similarity_keys.remove('gap')
        for key in similarity_keys:
            if not all(x in similarity[key] for x in similarity_keys) or \
                    len(similarity[key]) != len(similarity_keys):
                print("Invalid alignment similarity matrix")

        self.similarity = similarity
        if 'gap' not in self.similarity:
            self.similarity['gap'] = 0

        self.matrix = None
        self.seq1 = None
        self.seq2 = None


    def align(self, seq1, seq2):
        self.matrix = np.ndarray(shape=(len(seq1) + 1, len(seq2) + 1), \
                dtype=np.dtype([('score', np.int), ('trace', np.byte)]))
        seq1 = seq1.upper()
        seq2 = seq2.upper()

        # build matrix of scores/traceback information
        for i in xrange(len(seq1) + 1):
            self.matrix[i, 0] = (self.similarity['gap'] * i, Aligner._DEL)
        for j in xrange(len(seq2) + 1):
            self.matrix[0, j] = (self.similarity['gap'] * j, Aligner._INS)
        for i in xrange(1, len(seq1) + 1):
            for j in xrange(1, len(seq2) + 1):
                match = (self.matrix[i - 1, j - 1]['score'] + \
                            self.similarity[seq1[i - 1]][seq2[j - 1]], Aligner._MAT)
                delete = (self.matrix[i - 1, j]['score'] + \
                            self.similarity['gap'], Aligner._DEL)
                insert = (self.matrix[i, j - 1]['score'] + \
                            self.similarity['gap'], Aligner._INS)
                self.matrix[i, j] = max(delete, insert, match, 
                                  key=lambda x: x[0])
        self.matrix[0, 0] = (0, Aligner._END)

        # calculate alignment from the traceback
        i = len(seq1)
        j = len(seq2)
        traceback = list()
        while i > 0 or j > 0:
            if self.matrix[i, j]['trace'] == Aligner._MAT:
                if seq1[i - 1] == seq2[j - 1]:
                    traceback.append((i - 1, j - 1, "match", None))
                else:
                    traceback.append((i - 1, j - 1, "mismatch", None))
                i -= 1
                j -= 1
            elif self.matrix[i, j]['trace'] == Aligner._INS:
                traceback.append((i - 1, j - 1, "insertion", 1))
                j -= 1
            elif self.matrix[i, j]['trace'] == Aligner._DEL:
                traceback.append((i - 1, j - 1, "deletion", 1))
                i -= 1
            elif self.matrix[i, j]['trace'] == Aligner._END:
                pass
            else:
                # error
                pass
        traceback.reverse()

        # combine indels
        indel = None
        traceback_combined = list()
        for t in traceback:
            if t[2] == "insertion" or t[2] == "deletion":
                if indel is not None:
                    if t[2] == indel[2]:
                        indel[3] += t[3]
                    else:
                        print("Check gap penalty")
                else:
                    indel = list(t)
            else:
                if indel is not None:
                    traceback_combined.append(tuple(indel))
                    indel = None
                traceback_combined.append(t)
        if indel is not None:
            traceback_combined.append(tuple(indel))

        return traceback_combined


