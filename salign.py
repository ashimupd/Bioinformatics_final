'''
This is a dictionary of dictionary holding the values as it would be in BLOSUM62 substitution matrix.
This substiution matrix along with affine gap penalties is the basis of our scoring model
'''
blosum62 = {
'C':{'C':9, 'S':-1, 'T':-1, 'P':-3, 'A':0,  'G':-3, 'N':-3, 'D':-3, 'E':-4, 'Q':-3, 'H':-3, 'R':-3, 'K':-3, 'M':-1, 'I':-1, 'L':-1, 'V':-1, 'F':-2, 'Y':-2, 'W':-2},
'S':{'C':-1,'S':4,  'T':1,  'P':-1, 'A':1,  'G':0,  'N':1,  'D':0,  'E':0,  'Q':0,  'H':-1, 'R':-1, 'K':0,  'M':-1, 'I':-2, 'L':-2, 'V':-2, 'F':-2, 'Y':-2, 'W':-3},
'T':{'C':-1,'S':1,  'T':4,  'P':1,  'A':-1, 'G':1,  'N':0,  'D':1,  'E':0,  'Q':0,  'H':0,  'R':-1, 'K':0,  'M':-1, 'I':-2, 'L':-2, 'V':-2, 'F':-2, 'Y':-2, 'W':-3},
'P':{'C':-3,'S':-1, 'T':1,  'P':7,  'A':-1, 'G':-2, 'N':-1, 'D':-1, 'E':-1, 'Q':-1, 'H':-2, 'R':-2, 'K':-1, 'M':-2, 'I':-3, 'L':-3, 'V':-2, 'F':-4, 'Y':-3, 'W':-4},
'A':{'C':0, 'S':1,  'T':-1, 'P':-1, 'A':4,  'G':0,  'N':-1, 'D':-2, 'E':-1, 'Q':-1, 'H':-2, 'R':-1, 'K':-1, 'M':-1, 'I':-1, 'L':-1, 'V':-2, 'F':-2, 'Y':-2, 'W':-3},
'G':{'C':-3,'S':0,  'T':1,  'P':-2, 'A':0,  'G':6,  'N':-2, 'D':-1, 'E':-2, 'Q':-2, 'H':-2, 'R':-2, 'K':-2, 'M':-3, 'I':-4, 'L':-4, 'V':0,  'F':-3, 'Y':-3, 'W':-2},
'N':{'C':-3,'S':1,  'T':0,  'P':-2, 'A':-2, 'G':0,  'N':6,  'D':1,  'E':0,  'Q':0,  'H':-1, 'R':0,  'K':0,  'M':-2, 'I':-3, 'L':-3, 'V':-3, 'F':-3, 'Y':-2, 'W':-4},
'D':{'C':-3,'S':0,  'T':1,  'P':-1, 'A':-2, 'G':-1, 'N':1,  'D':6,  'E':2,  'Q':0,  'H':-1, 'R':-2, 'K':-1, 'M':-3, 'I':-3, 'L':-4, 'V':-3, 'F':-3, 'Y':-3, 'W':-4},
'E':{'C':-4,'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':2,  'E':5,  'Q':2,  'H':0,  'R':0,  'K':1,  'M':-2, 'I':-3, 'L':-3, 'V':-3, 'F':-3, 'Y':-2, 'W':-3},
'Q':{'C':-3,'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':0,  'E':2,  'Q':5,  'H':0,  'R':1,  'K':1,  'M':0,  'I':-3, 'L':-2, 'V':-2, 'F':-3, 'Y':-1, 'W':-2},
'H':{'C':-3,'S':-1, 'T':0,  'P':-2, 'A':-2, 'G':-2, 'N':1,  'D':1,  'E':0,  'Q':0,  'H':8,  'R':0,  'K':-1, 'M':-2, 'I':-3, 'L':-3, 'V':-2, 'F':-1, 'Y':2,  'W':-2},
'R':{'C':-3,'S':-1, 'T':-1, 'P':-2, 'A':-1, 'G':-2, 'N':0,  'D':-2, 'E':0,  'Q':1,  'H':0,  'R':5,  'K':2,  'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':-3, 'Y':-2, 'W':-3},
'K':{'C':-3,'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':-1, 'E':1,  'Q':1,  'H':-1, 'R':2,  'K':5,  'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':-3, 'Y':-2, 'W':-3},
'M':{'C':-1,'S':-1, 'T':-1, 'P':-2, 'A':-1, 'G':-3, 'N':-2, 'D':-3, 'E':-2, 'Q':0,  'H':-2, 'R':-1, 'K':-1, 'M':5,  'I':1,  'L':2,  'V':-2, 'F':0,  'Y':-1, 'W':-1},
'I':{'C':-1,'S':-2, 'T':-2, 'P':-3, 'A':-1, 'G':-4, 'N':-3, 'D':-3, 'E':-3, 'Q':-3, 'H':-3, 'R':-3, 'K':-3, 'M':1,  'I':4,  'L':2,  'V':1,  'F':0,  'Y':-1, 'W':-3},
'L':{'C':-1,'S':-2, 'T':-2, 'P':-3, 'A':-1, 'G':-4, 'N':-3, 'D':-4, 'E':-3, 'Q':-2, 'H':-3, 'R':-2, 'K':-2, 'M':2,  'I':2,  'L':4,  'V':3,  'F':0,  'Y':-1, 'W':-2},
'V':{'C':-1,'S':-2, 'T':-2, 'P':-2, 'A':0,  'G':-3, 'N':-3, 'D':-3, 'E':-2, 'Q':-2, 'H':-3, 'R':-3, 'K':-2, 'M':1,  'I':3,  'L':1,  'V':4,  'F':-1, 'Y':-1, 'W':-3},
'F':{'C':-2,'S':-2, 'T':-2, 'P':-4, 'A':-2, 'G':-3, 'N':-3, 'D':-3, 'E':-3, 'Q':-3, 'H':-1, 'R':-3, 'K':-3, 'M':0,  'I':0,  'L':0,  'V':-1, 'F':6,  'Y':3,  'W':1},
'Y':{'C':-2,'S':-2, 'T':-2, 'P':-3, 'A':-2, 'G':-3, 'N':-2, 'D':-3, 'E':-2, 'Q':-1, 'H':2,  'R':-2, 'K':-2, 'M':-1, 'I':-1, 'L':-1, 'V':-1, 'F':3,  'Y':7,  'W':2},
'W':{'C':-2,'S':-3, 'T':-3, 'P':-4, 'A':-3, 'G':-2, 'N':-4, 'D':-4, 'E':-3, 'Q':-2, 'H':-2, 'R':-3, 'K':-3, 'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':1,  'Y':2,  'W':11}
}

'''
Implement neighbor joining algorithm here 
'''    
def neighbor_joining(initial_distance_matrix):
    print("To be implemented")




'''
Jon Beck
A function to read a fasta file
Last modified: 26 August 2020
'''

'''
readfasta - a functionto reads a fasta file
Parameter: a filename that must be in fasta format.  The file is
assumed to have:
1. the first line a header line
2. arbitrary blank lines
3. every line (especially including the last) is terminated by a newline
   terminator (enter key)
4. no line has only spaces on it

Return: a list of lists. Each inner list will have two elements:
1. the sequence identifier, i.e., the characters between the leading ">"
   and the first space
2. the sequence, a single string of all the letters with no line terminators
'''
def readfasta(filename):
    result_list = []
    with open(filename, 'r') as infile:

        # process the first line, which must be a header line
        line = infile.readline()
        header_line = line.rstrip()
        label = header_line[1:].split(' ')[0]

        # initialize the sequence accumulator
        sequence = ''

        # process all the rest of the lines in the file
        for line in infile:
            line = line.rstrip()

            # ignore blank lines
            if line != '':

                # if it's a header line, finish the previous sequence
                # and start a new one
                if line[0] == '>':
                    result_list.append([label, sequence])

                    label = line[1:].split(' ')[0]
                    sequence = ''
            
                    # if we're here, we must be in letters of the sequence
                else:
                    sequence += line
            
    # we're done, so clean up, terminate the last sequence, and return
    result_list.append([label, sequence])
    return result_list

'''
Takes two protein sequences and calculates the score and traceback matrix
In order to calculate, the score matrix this algorithm uses BLOSUM62 substitution score matrix.
The Scoring model is called differential weighting of mismatch penalties.
g = -8 is used for gaps. This is the commonly used value for gaps in protein sequence alignment.
@param seq1: protein sequence 1 to be aligned  
@param seq2: protein sequence 2 to be aligned with sequence 1 
@param score: initial score of the given promoter sequence
@param g: gap introduction penalty
@return score: a matrix of alignment scores filled up using dp
@return trace: a matrix that will help in tracing the optimal alignment.
'''
def seq_align(seq1, seq2, g):


    m = len(seq1)
    n = len(seq2)
    # This is a score matrix
    score = [ [ 0 for i in range(n+1) ] for j in range(m+1) ]
    
    # Initialize the matrix according to the score model defined above
    for i in range(m+1):
        score[i][0] = g * i # The first column has the value of gap in increasing sequence
    
    for j in range(n+1):
        score[0][j] = g * j # The first row has the value of gap in increasing sequence
    
    
    # This is a traceback matrix
    # 1 inside the matrix means value came from diagonal.
    # 2 inside the matrix means value came from up.
    # 3 inside the matrix means that value came from left.

    trace = [ [ 0 for i in range(n+1) ] for j in range(m+1) ]
    
    # Initialize the matrix according to the traceback model defined above
    for i in range(m+1):
        trace[i][0] = 2 # value came from up, down the first colum
    
    for j in range(n+1):
        trace[0][j] = 3 # value came from left, along the first row
    
   
    # Filling the matrix using DP
    for i in range(0, m):
        for j in range(0,n):
            # value from diagonal means that it is either a match or unmatch.
            # Depending upon the types of unmatch, different score is assigned.
            # This is called differential weighting of unmatched pair.
            diag = score[i][j] + blosum62[seq1[i]] [seq2[j]] # value from the diagonal
            up = score[i][j+1] + g # value from up 
            left = score[i+1][j] + g # value from left
            score[i+1][j+1] = (max(diag, up, left)) # keep the maximum of three choices
            trace[i+1][j+1] = direction_finder(diag, up, left) # filling the traceback matrix with the correct value as defined above depending upon-  
                                                               # what was the maximum value among diag, up and left 
	
    #print(score)
    #print(trace)
    return(score, trace)
    
'''
This function traces back the score matrix with the help of traceback matrix and returns the optimal alignment
@param score: score matrix
@param trace: traceback matrix
@param seq1: the original protein sequence 1
@param seq2: the original protein sequence 2
@return seq1, seq2: optimal alignment of two sequences
'''
def trace_back(score, trace, seq1, seq2):
    
    align1 = "" # will store optimal alignment for seq1
    align2 = "" # will store optimal alignment for seq2
    
    i = len(seq1) # length of proetin sequence 1
    j = len(seq2) # length of protein sequence 2
    # Traceback starts from the bottom left column (last cell of the matrix)
    
    while i > 0 or j > 0:
        # this means respective characters (amino acids) are aligned from both strings, no gaps.
        if trace[i][j] == 1:
            align1 = seq1[i-1] + align1 # Saving optimal alignment in each iteration for protein seq 1
            align2 = seq2[j-1] + align2  # Saving optimal alignment in each iteration for protein seq 2.
            i = i - 1 # time for next character (amino acid) in the sequence 1
            j = j - 1 # time for next character (amino acid) in sequence 2
        # this introduces a gap in seq 2
        # this means a character in that poistion in seq 1 is aligned with a gap in seq 2
            
        elif trace[i][j] == 2:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i = i - 1
          
        # this introduces a gap in seq 1
        # this means a character in that poistion in seq 1 is aligned with a gap in seq 1
        else:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j = j - 1
            
    return(align1, align2)

'''
This is a helper function that finds the correct direction of where the value came from in the score matrix
In other words finds out which variable has the maximum value between 'diag', 'up', and 'left' 
@param diag: value from diagonal
@param up: value from up
@param left: value from left
@return 1: if 'diag' has the maximum value
@return 2: if 'up' has the maximum value
@return 3: if 'left' has the maximum value
''' 
def direction_finder(diag, up, left):
   
   if diag > up:
       if diag > left: return 1 # diag is greater than up and left both
       else: return 3           # diag is greater than up but less than left so left is the greatest
   else:
        if up > left: return 2  # up is greater than diag and left both
        else: return 3          # up is greater than diag, but less than left so left us the greatest.
 
'''
This function calculates the alignment score taking into consideration the affine gap penalty model
The alignment score is taken as an estimation of evolutionary distances.
This means an introduction of gap is penalized heavily than an extension of gap
@param g: gap initiation penalty 
@param r: gap extension penalty
@param sequences: a list of list containing two aligned sequences.
@return score: total alignment score calculated between two sequences which is taken as an estimate of evolutionary distance.
'''
def affine_gap_score_align(align1, align2, g, r):
    aa_compared = 0.0 # number of residue compared excluding gap
    score = 0.0 # stores the final alignment score
    # Intially, in gap tracker for both sequence will be false 
    seq1_in_gap = False
    seq2_in_gap = False
    # after alignment both will have equal lengths
    for i in range(len(align1)):
        if align1[i] == '-':
        # If gap is encountered, lets check if we are already inside gap or not
            if seq1_in_gap: # if we are already inside a gap, then only assign gap extension penalty
                            # The gap extension penalty is 'r' in our case
                score = score + r
            else: # if we are not inside gap, after this we will be.
                  # Because this line is reached only after we encounter a gap
                  # Set seq1_in_gap to true
                  # assign a gap initiation penalty 'g' because this is where gap was introduced  
                score = score + g
                seq1_in_gap = True
        
        # We repeat same thing for alignment 2
        # Remember the length of two alignments are always equal
        elif align2[i] == '-':
             # If gap is encountered, lets check if we are already inside gap or not
            if seq2_in_gap: # if we are already inside a gap, then only assign gap extension penalty
                            # The gap extension penalty is 'r' in our case
                score = score + r
            else: # if we are not inside gap, after this we will be 
                  # Because this line is reached only after we encounter a gap
                  # Set seq1_in_gap to true
                  # assign a gap initiation penalty 'g' because this is where gap was introduced  
                score = score + g
                seq2_in_gap = True
        
        # Now we handle the scoring for non-gap cases. We make the use of 'BlOSUM62' matrix for penalizing the mismatches-
        # in the final alignment scoring as well
        # The below line is reached only when we are not inside gap.
        # Therefore we set the gap tracker to False
        
        else:
            if seq1_in_gap:
                seq1_in_gap = False
            if seq1_in_gap:
                seq2_in_gap = False
            # Reaching this line means, there is aa in both sequences at this position and not gap.
            # We are counting the number of times non-gap comparison is made between the aligned sequences.
            aa_compared = aa_compared + 1
            
            # Each match between aa in aligned sequences contributes +1 to the score
            if align1[i] == align2[i]:
                score = score + 1       
    return round((1.0 - (score/aa_compared)), 4) # return the score upto 3 decimal precision

'''
This is a helper function that returns the evolutionary distance estimate between the sequences-
by taking two sequences
@param seq1: protein sequence 1
@param seq2: protein sequence 2
@param g: gap initiation penalty
@param r: gap extension penalty
@return p_dist: pairwise distance
'''
def pairwise_distance(seq1, seq2, g, r):
    score, trace = seq_align(seq1, seq2, g) # this will return a score matrix and a traceback matrix
    align1, align2 = trace_back(score, trace, seq1, seq2) # this will do the traceback and return two aligned sequences
    p_dist = affine_gap_score_align(align1, align2, g, r) # this will return the pairwise distance between two sequences
    return p_dist # pairwise distance between two sequences is returned

'''
Creates a matrix of pairwise distances between sequences after thay have been aligned
@param file: name of file to be read that contains sequences in fasta format
@param g: gap initiation penalty
@param r: gap extension penalty
@return distance_matrix: A matrix of pairwise distance between sequences
'''
def create_distance_matrix(file, g, r):
    # reads sequences from the file and stores in a list of lists 
    # The 1st element of each inner list is the name of sequence 
    # The 2nd element of each inner list is the sequence itself 
    sequences = readfasta(file) 
    
    distance_matrix = []  # declaring a distance matrix where we would store final result
    # selects one sequence from the "file" one a at a time
    for i in range (0, len(sequences)):
        # declaring a a temp list where we would store pairwise distances between-
        # each sequence with every other sequences
        # for example, in first iteration, 'temp_dist_list' will contain the distance between itself and all-
        # the other sequences.
        temp_dist_list = []

        # Also selects one sequence from the "file" one at a time
        for case in sequences:           
            # calling the pairwise_distance function 
            p_dist = pairwise_distance(case[1], sequences[i][1], g, r)
            # appending the distance to the temporary list
            temp_dist_list.append(p_dist)
        # appending the temp list to the final list
        distance_matrix.append(temp_dist_list)
    return distance_matrix  
    
'''
     In the order it is found in sequences_cytb.txt

     distance_matrix =          Gene1   Gene2   Gene3   .....
                         Gene1  0         9       6    .....
                         Gene2  9         0       7    .....
                         Gene3  6         7       0    .....                           
'''
'''
A helper function that will display the distance matrix in a formatted way for the eye
@param init_matrix: pairwised distance matrix
'''
def display_formatted(init_matrix):
    m = init_matrix
    for i in range(0, len(m)):
        for j in range(0, len(m[i])):
            print(str(m[i][j]) + "|"),
        print("\n")


def main():
    
    # Code takes some time to run with number alignments goining on
    
    gap = -4  # gap initiation penalty
    affine_gap = -2 # gap extension penalty
    ''' 
    seq1 = "KTNIRKTHPLAKIVNNSFIDLPTPSNISAWWNFGSLLGICLLLQIATGLFLAMHYTPDTTTAFSSVTHICRDVNYGWIIRYMHANGASMFFMCLFLHVGRGLYYGSYTFSETWNIGVILLLTVMATAFMGYVLPWGQMSFWGATVITNLLSAIPYIGSNLVEWIWGGFSVDKATLTRFFAFHFILPFIITTLAMIHLLFLHETGSNNPAGITSESDKIPFHPYYTIKDILGILFLIFMLMGLVLFTPDLLGDPDNYTPANPLNTPPHIKPEWYFLFAYAILRSIPNKLGGVLALVLSILILIIIPLLHTSKQRSMMFRPLSQCLFWFLVADLLVLTWIGGQPVEYPFIIIGQLASIFYFMILLVLMPTMNIIENNLLKW"
    seq2 = "KTNIRKTHPLAKIVNNSFIDLPTPSNISAWWNFGSLLGICLLLQIATGLFLAMHYTPDTTTAFSSVTHICRDVNYGWIIRYMHANGASMFFMCLFLHVGRGLYYGSYTFSETWNIGVILLLTVMATAFMGYVLPWGQMSFWGATVITNLLSAIPYIGSNLVEWIWGGFSVDKATLTRFFAFHFILPFIITTLAMIHLLFLHETGSNNPAGITSESDKIPFHPYYTIKDILGILFLIFMLMGLVLFTPDLLGDPDNYTPANPLNTPPHIKPEWYFLFAYAILRSIPNKLGGVLALVLSILILIIIPLLHTSKQRSMMFRPLSQCLFWFLVADLLVLTWIGGQPVEYPFIIIGQLASIFYFMILLVLMPTMNIIENNLLKW"
    s ,t = seq_align(seq1,seq2,gap)
    a1,a2= trace_back(s,t,seq1,seq2)
    print(a1)
    print(a2)
    d = pairwise_distance(seq1,seq2,gap,affine_gap)
    print(d)    
    '''
    # use this matrix as an input to Neighbor Joining
    dmatrix = create_distance_matrix("sequences_cytb.txt", gap, affine_gap)
    # distances are rounded to 4th percision. If more precision is required can be changed in affine_gap_score_align(...) function
    
    display_formatted(dmatrix) # display pairwise distance matrix in a way that could be seen by eyes.
   
   
if __name__ == '__main__':
    main()






		
       
