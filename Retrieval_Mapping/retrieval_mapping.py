#!/usr/bin/env python
"""
Algorithms in Bioinformatics
Author: Diego Montiel
Student number: 880505580110
Implementation of the SSAHA algorithm

Hints:
- write a function to generate a hash table of k-mers from a set of sequences
- write a function to return a sorted list of hits in the hash table 
for a query sequence
- write a function to find the best hit
- write a fasta parse_fasta_file to read in the Arabidopsis data

"""
# import statements
from collections import Counter, defaultdict
from sys import argv
import time

# implement your functions here

def get_hash_table(seqs, kmer):
    """
    Construct a Hash table from the reference genome
    Loop the sequences and split in the number of kmer
    
    Output: 
        A hashtable dictionary of kmers (W) and as value tuple(s) 
    of integer values; Number of sequence and the positions 
    of W in the sequence

    Input: 
        seqs - List of N sequences
        kmer - integer of the length partition 
    
    """
    hash_table = {}
    #get all the index where k does not pass over the text
    for index, seq in enumerate(seqs):
        #start from zero, till len of seq, every kmer step
        for k in range(0,len(seq) - kmer + 1,kmer):
            #get  the kmer and non overlapping k tuples
            kmer_seq   = seq[k:k+kmer]
            if kmer_seq in hash_table:
                hash_table[kmer_seq] += [(index+1,k+1)]
            else:
                hash_table[kmer_seq] = [(index+1,k+1)]
        
    return hash_table

def get_list_m(query,hash_table,kmer):

    """
    Construct a master list M of accumulated hits
    
    Output: 
        A list of integer tuples; 
            1) index - number of the sequence
            2) Shift 
            3) Offset - position of the kmer in the sequence 
      
    Input: 
        seqs - List of N sequences
        kmer - integer of the length partition 
    """

    h = []
    for t in range(len(query) - kmer +1):
        if query[t:t+kmer] in hash_table:
            positions = hash_table[query[t:t+kmer]]
            for index in positions:
                h.append((index[0],index[1]-t,index[1]))
    m = sorted(h)

    return m

def get_matches(m):

    """
    Build a list of tuples of all the matches.
    Condition: at least one exact match in:
        index  = index  + 1 
        shift  = shift  + 1
    
    Output: 
        List of dictionary of matches from table M; 
            1) index - int value number of the sequence
            2) Shift - list of all the shifts 
            3) Offset - list of all the position of kmer position 
    Input: 
        M - List of tuples of hits of query and hashtable
    """
    #Iterate M to get all matches
    #in a sequence
    matches = []
    offset = 0
    counter = 0
    #iterate with the length of m
    for length in range(len(m)-1):
        if m[length][0] == m[length+1][0] and m[length][1] == m[length+1][1]:
            if counter < 1:
                # if m[length][index] is greater than offset it is a new match
                # so we append it to the list ot create new element
                if m[length][0] > offset:
                    matches.append({"index":int,\
                                       "shift":[], \
                                       "offset":[],\
                                       "counter":int})
                #Get always the index of the last append in the list
                index = len(matches)-1
                matches[index]["index"]  =  m[length][0]
                matches[index]["shift"] += [m[length][1]]
                matches[index]["offset"]+= [m[length][2]]

                matches[index]["index"]  =  m[length+1][0]
                matches[index]["shift"] += [m[length+1][1]]
                matches[index]["offset"]+= [m[length+1][2]]
                counter  = 1
                offset   = m[length+1][0]
            else:
                index = len(matches)-1  

                matches[index]["index"]  =  m[length+1][0]
                matches[index]["shift"] += [m[length+1][1]]
                matches[index]["offset"]+= [m[length+1][2]]
        else:
            counter = 0
    return matches

def get_longest_match(matches):

    """
    Get the longest match between the query and the database 
    (sequences).
    Condition: Look for the highest repeat in the sequence index
    and the shift :
    
    Output: 
        List of dictionary(ies) of the longest matches from a query 
            1) index   - int value number of the sequence
            2) Shift   - list of all the shifts 
            3) Offset  - list of all the position of kmer position 
            4) Counter - int of how many times is match repeated
    Input: 
        Matches - List of dictionary(ies) of matches of a query and a hashtable
    """
    result_max     = []
    longest_match  = matches
    length_matches = len(longest_match)
    for init in range(0,length_matches):
        #Get the maximal repeats from the shift
        d = defaultdict(int)
        for j in longest_match[init]["shift"]:
            d[j] += 1 
            result = max(d.iteritems(), key=lambda x: x[1])
   
         #Get the longest match per sequence
        result_max.append(result)    
        max_long_match = result_max[init][0]
        for value in range(len(longest_match[init]["shift"])-1,-1,-1):
            
            if longest_match[init]["shift"][value] != max_long_match:
                del longest_match[init]["shift"][value]
                del longest_match[init]["offset"][value]
                longest_match[init]["counter"] = \
                len(longest_match[init]["shift"])
    #Get the longest match sequence 
    if len(longest_match) > 1:
        for index in range(len( longest_match)-1,-1,-1):
            if longest_match[index]["counter"] <  \
               longest_match[index-1]["counter"]:
                del longest_match[index]
            elif longest_match[index]["counter"] > \
                 longest_match[index-1]["counter"]:
                del longest_match[index-1]

    return longest_match

def print_alignment(query, long_match, hash_table,kmer):

    """
    function that print the alignment of the query and the subject
    of the match. (two sequences) :
    
    Output: 
       Print the query, alignment and sequence
    Input: 
        Query - query sequence string
        Long_match - list of dictionaries with the longest matches
        hash_table - hash table from the database (sequences)
        kmer - size of the kmer sequence partition
    """

    tuple_match = []
    query_match = []
    sequence    = []
    align       = []
    
    #Convert the values from the long match list in a tuple
    for l in range(len(long_match)):
        for o in range(len(long_match[l]["offset"])): 
            tuple_match.append((long_match[l]["index"],long_match[l]["offset"][o]))
            seq = long_match[l]["index"]
    
    #Get the indexes from a hashtable
    sequence_pos = []
    for key,value in hash_table.items():
        for v in value:
            if v[0] == seq:
                sequence_pos.append(v[1])
                
    sequence_pos = sorted(sequence_pos) 

    #Get the alignments
    for index in sequence_pos:
        t = (seq, index)
        for key,value in hash_table.items():
            for v in value:
                if v == t:
                    sequence.append( key )
                    if v in tuple_match:  
                        query_match.append( key )
                        for j in range(kmer):
                            align.append("|")                            
                    else:
                        for j in range(kmer):
                            query_match.append( " " )                    
                            align.append(" ")
                
    print "Alignment:"
    print ''.join(query_match)
    print ''.join(align)
    print ''.join(sequence)

def parse_fasta_file(filename):
    """-Function: parseFastaFile
    Description: 
        Function that parse a FASTA file
    Input: 
        A fasta file
    Output: 
        Returns a list of dictionaries:
           Including the label, sequence and the lenght sequence
    """
    fasta_file      = open(filename,'r')
    tmp_seq         = []
    sequence        = []
    name            = []
    length          = []
    index           = 0
    for line in fasta_file:
        if line[0] == '>':
            name += [line[1:].rstrip()]
            if tmp_seq != []:
                tmp = ''.join(tmp_seq)
                sequence.append(tmp)
                length.append(len(sequence[index]))
                tmp_seq = []
                index  +=1
        else:
            line = line.rstrip()
            tmp_seq.append(line)
    tmp = ''.join(tmp_seq)
    sequence.append(tmp)
    tmp_seq    = []
    fasta_file.close()
    return name, sequence, length

if __name__ == "__main__":

    # the code below should produce the results necessary to answer the questions
    # in other words, if we run your code, we should see the data that you used
    # to answer the questions

    #####################
    # Initial Settings #
    # The script in order to work needs three parameters
    # 1) Fasta file of the sequence genome
    # 2) Fasta file of the sequence query(ies)
    # 3) Integer value of the kmer
    #seq_fasta   = argv[1]
    #query_fasta = argv[2]
    #kmer        = argv[3]
    #########################

    start_time  = time.time()

    ###########################
    # Test set
    ###########################
    kmer         = 2    
    query       = "TGCAACAT"
    s1          = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
    s2          = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
    s3          = 'GGATCCCCTGTCCTCTCTGTCACATA'
    seqs        = [s1,s2,s3]
    hash_table  = get_hash_table(seqs,kmer)
    #hash_table  = builderHT(seqs,kmer)
    m           = get_list_m(query,hash_table,kmer)
    matches     = get_matches(m)

    print "##############"
    print "# Question 1:#"
    print "##############"
    print "Hash table"
    keylist = hash_table.keys()
    keylist.sort()
    for key in keylist:
       print "%s: %s" % (key, hash_table[key])
    print " \n" 

    print "##############"
    print "# Question 2:#"
    print "##############"
    print "Number of hits: ",len(m)
    print "Matches: First and Last position "
    for index in range(len(matches)):
        print matches[index]["index"],min(matches[0]["shift"]),\
              min(matches[0]["offset"])
        print matches[index]["index"],max(matches[0]["shift"]),\
              max(matches[index]["offset"])
    long_match  = get_longest_match(matches)
    print " \n" 
    print "##############"
    print "# Question 3:#"
    print "##############"
    print_alignment(query, long_match, hash_table,kmer)   
    print "\n"
    print "Loading please wait..."

    ################################
    # Arabidopsis Genome
    ###############################

    kmer        = 15
    seq_fasta   = "TAIR10.fasta"
    query_fasta = "athal_query.fasta"

    names, seqs, length = parse_fasta_file(seq_fasta)
    q_name, query, _    = parse_fasta_file(query_fasta)
    hash_table          = get_hash_table(seqs,kmer)
    
    print "##############"
    print "# Question 4 "
    print "##############"
    print "Sequences:"
    print names
    print "Total: ", len(names)
    print "\n"
    print "Lengths:"
    print length
    print "Total: ", sum(length)

    print "\n"
    print "################"
    print "# Question 5: #"
    print "################"
    print "Total of 15-mers in the hash table:", len(hash_table.keys())
    print "\n"

    for i in range(len(query)):  
        m           = get_list_m(query[i],hash_table,kmer)
        matches     = get_matches(m)
        long_match  = get_longest_match(matches)

        print "################################"
        print "# Question 6 - Sequence: #", q_name[i]
        print "################################"
        print "Maximum of hit(s)"
        count = 0
        for index in range(len(long_match)):
            for i in range(len(long_match[index]["shift"])):
                print long_match[index]["index"],"", \
                      long_match[index]["shift"][i]\
                ,"",long_match[index]["offset"][i]
                count+=1
        print "Number of hits: ",count 
   
    print("--- Script has finished in %s seconds ---" \
          % (time.time() - start_time))
    
