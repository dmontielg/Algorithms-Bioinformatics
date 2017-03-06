#!/usr/bin/env python

"""
Author: Diego Montiel
Student nr: 880505580110
Script to: Network inference and enhancement
Boolean Networks REVEAL Algorithm
"""

from __future__ import division
from itertools import combinations
from copy import deepcopy
import math
import time

eps = 1e-100

def log2(x):
    """
    Function that returns the log2 of a int or float value
    
    Input: Int or float value
    """
    return math.log(x+eps)/math.log(2)

def entropy_single(x):
    """
    Function that returns the entropy of a single variable
    
    Input: List of boolean values from a probability
    """
    dictionary = {0:0,1:0}
    pi = 0
    pj = 0
    h = 0    
    for id,value in enumerate(x):
        dictionary[value] += 1
    pi = dictionary[0]/len(x)
    pj = dictionary[1]/len(x)
    h = -pi*log2(pi)-pj*log2(pj)

    return round(h,3)

def entropy_multiple(args):
    """
    Function that returns the entropy of H(X') and H(X,Y)
    
    Input: Lists of boolean probabilities
    """
    dictionary = {}
    probs   = {}
    xyz     = ''
    pij     = 0
    h       = 0 

    for idx,value in enumerate(args):
            dictionary[idx] = value
    for value in range(len(dictionary[0])):
        xyz = ''
        for idx in range(len(dictionary.keys())):
            #Concatenate each of the column of each list of bool
            #to get a certain combination, e.g. 001, 002
            xyz += str(dictionary[idx][value])
        if xyz in probs:
            #If already exist that combination in the list just add 1
            probs[xyz] += 1
            #if not create the combination in the list and add 1
        else:
            probs[xyz] = 1
            
    for value in probs.values():    
        pij = value/ sum(probs.values())
        h += -pij*log2(pij)
    return round(h,3)

def entropy(*args):
    """Return the entropy H(x,y) of x, or the joint entropy of N number of variables

    Input: 1 list of boolean values or a tuple of lists of boolean values
    """
    
    if len(args) == 1:
        for value in args:
            tup = value
        #if is a list means only is a single list
        if type(tup) is list:        
            return entropy_single(args[0])
        #if is a tuple means is a multiple list 
        else:
            return entropy_multiple(tup) 
    else:
       return entropy_multiple(args)   

def mutual_information_equation(args):
    """Return the M value of the formal equation of the mutual information

    Input: Tuple of lists of booleans
    """

    h_x    = args[0]
    h_yzw  = args[1:]
    h_xyzw = args
    m = entropy(h_x) + entropy(h_yzw) - entropy(h_xyzw)          

    return m

def mutual_information(*args):
    """Return the mutual information M(x,y) between x and y,
       or X against N number of information.

    Input: Value of the mutual information e.g of X against and Y
    """
    if len(args) < 2:
        print "Too few parameters, you need at least two!"
        return False
    else:
        tup = None
        if type(args[1][0]) is list:
            for idx,value in enumerate(args[1:]):
                tup = (args[0],)
                for v in value:
                    tup += (v,)       
            args = tup
            return mutual_information_equation(args)
        else:
            return mutual_information_equation(args)

def reveal(input,output,kmax=3):
    """Returns network, a dictionary containing for every target node
       a list of most likely source nodes. kmax determines how large
       the subsets explored are.

       Input:  list of booleans from input and output
    """
    bool_network = {}
    key     = sorted(input.keys())
    #Step 1) Pairwise comparison
    validation = []  
    for idx,value in enumerate(key):
        element = key[idx][0]
        h_x_prime = entropy(output[element])
        #Loop all of the x values (keys)
        for h in key:
            h_x_prime_x = entropy(output[element],input[h])
            #if H(X',Y) == H(X')
            if abs(h_x_prime - h_x_prime_x) < 0.001:
                bool_network[element] = h
                index = key.index(element)
                validation.append(key[index])
    if kmax == 1:
        return bool_network
    if len(key) == len(validation):
        return bool_network
    else:
        #Step i) check all the combinations given a kmax
        #starting from 2 and so on until reach the valueo of kmax
        k = 2
        while k <= kmax:
            #Gets all the combinations from the keys
            combination_k = list(combinations(sorted(input.keys()),k))
            #Start from each of the value hx
            for idx,value in enumerate(key):
                hx = entropy(output[value])
                #Check against all the possible combinations
                if value in validation:
                    pass
                else:
                    for comb in combination_k:
                        comb = list(comb)
                        combination = [input[x] for x in comb]
                        mutual_xy = mutual_information(output[value],\
                                                combination)
                        #if M(X[X,Y]) == H(X')
                        if abs(hx - mutual_xy) < 0.001:
                            bool_network[value] = comb
                            index = key.index(value)
                            validation.append(key[index])
                            break
            k += 1
        return bool_network

def read_tsv_file(tsv):
    """Reads a timeseries dataset and returns the dictionaries
       inputs (time 0,...,T-1) and outputs (time 1,...,T).
    """
    
    header = tsv.readline().rstrip()
    nodes  = header.split("\t")
    
    inputs  = {}; outputs = {}
    for i in range(1,len(nodes)):
      inputs[nodes[i]] = []; outputs[nodes[i]] = []

    line = tsv.readline().rstrip()
    
    while len(line) > 0:
        vals = line.split("\t")
        for i in range(1,len(nodes)):
            inputs[nodes[i]].append(float(vals[i]))
        line = tsv.readline().rstrip()

    outputs = deepcopy(inputs)
    n = len(inputs[nodes[1]])
    for i in range(1,len(nodes)):
    
        del inputs[nodes[i]][n-1];
        del outputs[nodes[i]][0]
          
    return inputs, outputs
                               
def check_table(input,output):
    """Prints all entropy and mutual information values
       in Figure 5 of Liang et al, PSB 1998 for checking
       your implementation.
    """
    print "H(A), H(B), H(C) =",       entropy(input['A']), entropy(input['B']), entropy(input['C'])
    print "H(A,B), H(B,C), H(A,C) =", entropy(input['A'],input['B']), entropy(input['B'],input['C']), entropy(input['A'],input['C'])
    print "H(A,B,C) =",               entropy(input['A'],input['B'],input['C'])

    print
    print "H(A') =", entropy(output['A'])
    print "H(A',A), H(A',B), H(A',C) =", \
        entropy(output['A'],input['A']), \
        entropy(output['A'],input['B']), \
        entropy(output['A'],input['C'])
    print "M(A',A), M(A',B), M(A',C) =", \
        mutual_information(output['A'],input['A']), \
        mutual_information(output['A'],input['B']), \
        mutual_information(output['A'],input['C'])
    print
    print "H(B') =", entropy(output['B'])
    print "H(B',A), H(B',B), H(B',C) =", \
        entropy(output['B'],input['A']), \
        entropy(output['B'],input['B']), \
        entropy(output['B'],input['C'])
    print "H(B',[A,B]), H(B',[B,C]), H(B',[A,C]) =", \
        entropy(output['B'],input['A'],input['B']), \
        entropy(output['B'],input['B'],input['C']), \
        entropy(output['B'],input['A'],input['C'])
    print "M(B',A), M(B',B), M(B',C) =", \
        mutual_information(output['B'],input['A']), \
        mutual_information(output['B'],input['B']), \
        mutual_information(output['B'],input['C'])
    print "M(B',[A,B]), M(B',[B,C]), M(B',[A,C]) =", \
        mutual_information(output['B'],input['A'],input['B']), \
        mutual_information(output['B'],input['B'],input['C']), \
        mutual_information(output['B'],input['A'],input['C'])

    print
    print "H(C') =", entropy(output['C'])
    print "H(C',A), H(C',B), H(C',C) =", \
        entropy(output['C'],input['A']), \
        entropy(output['C'],input['B']), \
        entropy(output['C'],input['C'])
    print "H(C',[A,B]), H(C',[B,C]), H(C',[A,C]) =", \
        entropy(output['C'],input['A'],input['B']), \
        entropy(output['C'],input['B'],input['C']), \
        entropy(output['C'],input['A'],input['C'])
    print "H(C',[A,B,C]) =", entropy(output['C'],input['A'],input['B'],input['C'])
    print "M(C',A), M(C',B), M(C',C) =", \
        mutual_information(output['C'],input['A']), \
        mutual_information(output['C'],input['B']), \
        mutual_information(output['C'],input['C'])
    print "M(C',[A,B]), M(C',[B,C]), M(C',[A,C]) =", \
        mutual_information(output['C'],input['A'],input['B']), \
        mutual_information(output['C'],input['B'],input['C']), \
        mutual_information(output['C'],input['A'],input['C'])
    print "M(C',[A,B,C]) =", mutual_information(output['C'],input['A'],input['B'],input['C'])
                                                                                                                                                                                                                      
if __name__ == "__main__":

    start_time  = time.time()
    # Data used in Liang et al, PSB 1998, in Fig. 1 and on pp. 23:
    
    input = {'A':[0,0,0,0,1,1,1,1],'B':[0,0,1,1,0,0,1,1],'C':[0,1,0,1,0,1,0,1]}
    output = {'A':[0,0,1,1,0,0,1,1],'B':[0,1,0,1,1,1,1,1],'C':[0,0,0,1,0,1,1,1]}

    print "###############"
    print "# Quesrion 1  #"
    print "###############"
    
    check_table(input,output)
    network = reveal(input,output)   
    print 
    for node in sorted(network.keys()):
        print "The expression of node", node, "is best explained by nodes", network[node]
    # Data simulated from a real network:
    print "Loading please wait..."
    input, output = read_tsv_file(open('yeast_bin.tsv'))
    network = reveal(input,output,1)
    
    print "###############"
    print "# Quesrion 3  #"
    print "###############"
    print
    for node in network.keys():
        print "The expression of node", node, "is best explained by nodes", network[node]
        
    print "###############"
    print "# Quesrion 5  #"
    print "###############"
    print 
    in_dics = {}
    out_dics = {}
    print "100 timepoints"
    for i in input:
        in_dics[i] = input[i][8999:9999]
        out_dics[i] = output[i][8999:9999]                    
    network = reveal(in_dics,out_dics,3)
    print
    for node in network.keys():
        print "The expression of node", node, "is best explained by nodes", network[node]
    print 
    print "1000 timepoints"
    for i in input:
        in_dics[i] = input[i][0:1000]
        out_dics[i] = output[i][0:1000]                    
    network = reveal(in_dics,out_dics,3)
    print
    for node in network.keys():
        print "The expression of node", node, "is best explained by nodes", network[node]
    print("--- Script has finished in %s seconds ---" \
          % (time.time() - start_time))
