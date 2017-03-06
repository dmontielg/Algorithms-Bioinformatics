#!/usr/bin/env python

"""
Algorithms in Bioinformatics
STUDENT NUMBER: 880505580110
Author:         Diego Montiel
"""
import time
import operator

#from sys import argv
#import random

##################################
# BLOSUM62 MATRIX
##################################
blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""

################################
# FUNCTIONS AND METHODS
################################

def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 23:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            blosum_matrix.append(map(int,parts[1:]))
    return order, blosum_matrix

def score(res1, res2):
    """
    Return similarity score from BLOSUM62 matrix for two residues
    
    res1: string, amino acid
    res2: string, amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1]
    lookup2 = BLOSUM62_ORDER[res2]
    return BLOSUM62_MATRIX[lookup1][lookup2]

def getAlphabet():
    """
    Generate ALPHABET (AMINO ACIDS)
    from BLOSUM62_ORDER
    """
    sorted_key = sorted(BLOSUM62_ORDER.items(), \
                        key=operator.itemgetter(1))
    ALPHABET = []
    for index, _ in sorted_key:
        ALPHABET.append(index)
    return ALPHABET

def arrayAsMatrix(matrix):
    """
    Returns a matrix from any list or array
    """
    return ('\n'.join([''.join(['{:3}'.format(item) for item in row])\
    for row in matrix]))

def initialMatrix(seq1,seq2):
    
    """
    Return an initial matrix where first row and first column is
    filled by the gap penalty
    
    seq1: string, protein sequence or amino acids
    seq2: string, protein sequence or amino acids    
    """
    initial_matrix = []
    #Construct a matrix of zeros with the length of both sequences
    for _ in range(len(seq1)+1):
        initial_matrix.append([0] * (len(seq2)+1))
        
    for row in range(1, len(seq1)+1):
        initial_matrix[row][0] = initial_matrix[row-1][0] + GAP_END
    
    for col in range(1, len(seq2)+1):
        initial_matrix[0][col] = initial_matrix[0][col-1] + GAP_END
    return initial_matrix

def globalAlignment(seq1,seq2):

    """
    #Return
    matrix: The optimal global alignment of two sequences, seq1
    and seq2   
    align_lenght: lenght size of the complete alignment between
    two sequences
    
    #Input
    seq1: string, protein sequence or amino acids
    seq2: string, protein sequence or amino acids 
    """

    nrow = len(seq1)
    ncol = len(seq2)

    #Create an initial matrix
    initial_matrix  = initialMatrix(seq1,seq2)
    matrix          = initialMatrix(seq1,seq2)
    #Fill the matrix with comparing each of the aligns
    for row in range(1, nrow):
        for col in range(1, ncol):

            distVer  = matrix[row-1][col]    + GAP_OPEN
            distHor  = matrix[row][col-1]    + GAP_OPEN
            distDiag = matrix[row-1][col-1]  + \
                BLOSUM62_MATRIX[ALPHABET.index(seq1[row-1])]\
                [ALPHABET.index(seq2[col-1])]
            #Stores the Maximal value from the three edit distances, 
            #horizontal, vertical and diagonal in the matrix.
            matrix[row][col] = max(distHor, distVer, distDiag)

    #Asigns the end gap to the last row
    for row in range(1, nrow):
        col = ncol
        matrix[row][col] = matrix[row-1][0]  + GAP_END    

    #Asigns the end gap to the last column
    for col in range(1, ncol):
        row = nrow
        matrix[row][col] = matrix[0][col-1]  + GAP_END    


    #Calculate the last row values according to the end gap
    #and open gap
    for row in range(1, nrow+1):
        col = ncol
        distHor  = matrix[row-1][col]    + GAP_END
        distVer  = matrix[row][col-1]    + GAP_OPEN
        distDiag = matrix[row-1][col-1]  + \
                BLOSUM62_MATRIX[ALPHABET.index(seq1[row-1])]\
                [ALPHABET.index(seq2[col-1])]
        matrix[row][col] = max(distHor, distVer, distDiag)

    #Calculate the last column values according to the end gap
    #and open gap
    for col in range(1, ncol+1):
        row = nrow
        distHor  = matrix[row-1][col]    + GAP_OPEN
        distVer  = matrix[row][col-1]    + GAP_END
        distDiag = matrix[row-1][col-1]  + \
                BLOSUM62_MATRIX[ALPHABET.index(seq1[row-1])]\
                [ALPHABET.index(seq2[col-1])]
        matrix[row][col] = max(distHor, distVer, distDiag)
        
    #Computing the final score of the matrix
    row = nrow
    col = ncol
    distVer  = matrix[row-1][col]    + GAP_END
    distHor  = matrix[row][col-1]    + GAP_END
    distDiag = matrix[row-1][col-1]  + \
                BLOSUM62_MATRIX[ALPHABET.index(seq1[row-1])]\
                [ALPHABET.index(seq2[col-1])]
    #calculate the maximal score value from the matrix
    matrix[row][col] = max(distHor, distVer, distDiag)
    #Call the function traceback to get the path of the alignment
    _, alignments,gaps = traceback(seq1,seq2,matrix)
    
    return matrix[-1][-1], alignments,initial_matrix, matrix,gaps

def traceback(seq1,seq2,matrix):

    """
    Returns 
    Array: Traceback of the global alignment with symbols
        ^ goes up, < goes left  and / goes in diagonal
    Int: lenght of the alignment
    
    Inputs: 
    seq1: amino acid sequences
    seq2: amino acid sequences    
    matrix: global alignment matrix
    """
    #Print score, which is last element in the matrix.
    traceback   = []
    alignment_1 = []
    alignment_2 = []
    alignments =  []
    row         = len(seq1)
    col         = len(seq2)
        
    while col > 0 or row > 0:
        
        if col == 0:
            #if index is in the first column
            traceback.append("^")
            row = row  -1
            alignment_1.append(seq1[row])
            alignment_2.append('-')
            
        elif row == 0:
            #if row is in the first row
            traceback.append("<")
            col = col -1
            alignment_1.append('-')
            alignment_2.append(seq2[col])
            
        #These two functions check if the value above is the right
        #path two follow by substracting and compare the result
        #with the end gap and open gap
        elif matrix[row][col] - matrix[row-1][col] == GAP_END:
            traceback.append("^")
            row = row -1
            alignment_1.append(seq1[row])
            alignment_2.append('-')

        elif matrix[row][col] - matrix[row-1][col] == GAP_OPEN:
            traceback.append("^")
            row = row -1
            alignment_1.append(seq1[row])
            alignment_2.append('-')

        elif matrix[row][col] - matrix[row][col-1] == GAP_OPEN:
            traceback.append("<")
            col = col -1
            alignment_1.append('-')
            alignment_2.append(seq2[col])

        elif matrix[row][col] - matrix[row][col-1] == GAP_END:
            traceback.append("<")
            col = col -1
            alignment_1.append('-')
            alignment_2.append(seq2[col])
            
        else:
            traceback.append('/')
            row = row -1
            col = col -1
            alignment_1.append(seq1[row])
            alignment_2.append(seq2[col])

    alignment_1.reverse()
    alignment_2.reverse()
    alignments.append(''.join(alignment_1))
    alignments.append(''.join(alignment_2))
    gaps = ''.join(alignment_1).count('-') +\
           ''.join(alignment_2).count('-')
    #print traceback
    return traceback, alignments,gaps

def calculateIdentity(alignments):

    """
    Return the identity percentage between two align sequences
    alignments: list of two alignment sequences
    """
    iterator = 0
    for index in range(len(alignments[0])):
        if alignments[0][index] == alignments[1][index]:
            iterator += 1
    return round(float(iterator)/len(alignments[0])*100,2)

##############################
# GLOBAL VARIABLES
##############################

BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()
ALPHABET = getAlphabet()
#change the value for open gap and end gap
GAP_OPEN = -4
GAP_END  = -4

#############################
# MAIN METHOD
#############################
if __name__=="__main__":

    """
    Set the alignments to the function globalAlignment(e.g. se1,seq2W)
    and and run the script.
    """
    start_time = time.time()
    seq1 = "THISLINE"
    seq2 = "ISALIGNED"
    seq3 = "MGLLCSRSRHHTEDTDENTQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQIKLLFQTGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDSAKYMLSSESIAIGEKLSEIGGRLDYPRLTKDIAEGIETLWKDPAIQETCARGNELQVPDCTKYLMENLKRLSDINYIPTKEDVLYARVRTTGVVEIQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHLFEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPCFEKTSFMLFLNKFDIFEKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRVFKIYRTTALDQKLVKKTFKLVDETLRRRNLLEA"
    seq4 = "MGSSCSRSHSLSEAETTKNAKSADIDRRILQETKAEQHIHKLLLLGAGESGKSTIFKQIKLLFQTGFDEAELRSYTSVIHANVYQTIKILYEGAKELSQVESDSSKYVISPDNQEIGEKLSDIDGRLDYPLLNKELVLDVKRLWQDPAIQETYLRGSILQLPDCAQYFMENLDRLAEAGYVPTKEDVLYARVRTNGVVQIQFSPVGENKRGGEVYRLYDVGGQRNERRKWIHLFEGVNAVIFCAAISEYDQMLFEDETKNRMMETKELFDWVLKQRCFEKTSFILFLNKFDIFEKKIQKVPLSVCEWFKDYQPIAPGKQEVEHAYEFVKKKFEELYFQSSKPDRVDRVFKIYRTTALDQKLVKKTFKLIDESMRRSREGT"
    #   Change which sequence want to alignment in globalAlignment.
    #   Hint. comment the print of initial and filled matrix
    #   if you choose long sequences
    alg_score, alignments,initial_matrix,matrix,gaps =\
               globalAlignment(seq1,seq2)

    identity = calculateIdentity(alignments)
    print 'Initial matrix: '
    print arrayAsMatrix(initial_matrix)
    print 'Filled matrix: '
    print arrayAsMatrix(matrix)
    print 'GAP OPEN: '+str(GAP_OPEN)
    print 'GAP END: ' +str(GAP_END)
    print 'Alignment lenght: '  +str(len(alignments[0]))
    print 'Identity: '+str(identity)+'%'
    print 'Gaps: '    +str(gaps)
    print 'Score: '   +str(alg_score)
    print 'Sequence 1: '+alignments[0]
    print 'Sequence 1: '+alignments[1]
    #Calculate the time to finish to run
    print("--- %s seconds ---" % (time.time() - start_time))
