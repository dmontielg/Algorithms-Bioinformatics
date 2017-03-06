#!/usr/bin/env python

"""
Assignment: 5
Author: Diego Montiel
Student nr: 880505580110
Script to: Graphs II network mining and module identification
"""

from __future__ import division
from numpy import sqrt, array, zeros, mean, delete, vstack, insert
from collections import defaultdict
from heapq import *


def distance(edges):
    """
    Function that returns a matrix with the node distances
    
    Input: List of tuples of the edges between nodes
    """
    dim = []
    for i in edges:
        dim.append(i[0])
    dim = set(((dim)))
    dim = max(dim)
    dim += 1
    #create a matrix of zeros with the dimension according the nodes
    distance_matrix = zeros((dim,dim))

    for col in range(0,dim):
        for row in range(0,dim):
            distance_matrix[col][row] = shortest_path(edges,col,row)[0]
    
    return distance_matrix

def association(distance_matrix):
    """
    Function that returns a matrix with the association of each position
    
    Input: matrix with the node distances of the graph
    """

    if type(distance_matrix) is list:
        dim = len(distance_matrix) 
    else:
        dim = distance_matrix.shape[0]
    association_matrix = zeros((dim,dim))
    for col in range(0,dim):
        for row in range(0,dim):
            d = distance_matrix[col][row]
            if d == 0:
                d = 1 
            association_matrix[col][row] = round(1/(d**2),4)
    return association_matrix
    
def correlation(x,y):
    """
    Function that returns pearson correlation 
    coefficient of two vectors
    
    Input x: Vector of numbers
    Input y: Vector of numbers
    """

    n   = len(x)
    xy  = 0
    x_x = 0
    y_y = 0
    sum_X = sum(x)
    sum_y = sum(y)

    for index in range(len(x)):
        
        xy  += x[index] * y[index]
        x_x += x[index] * x[index]
        y_y += y[index] * y[index]

    r = ((n*xy) - (sum_X*sum_y)) \
    / (sqrt((n*x_x) - (sum_X**2)) * sqrt((n*y_y) - (sum_y**2)))  
    
    return round(r,4)

def dissimilarity_matrix_from_correlation(association_matrix):
    """
    Function that returns a  dissimilarity matrix for each position
    with the pearson correlation coefficient
    
    Input: matrix with the associations
    """

    dim = association_matrix.shape[0]
    dissimilarity_matrix = zeros((dim,dim))
    for col in range(0,dim):
        for row in range(0,dim):
            x = list(association_matrix[col,:])
            y = list(association_matrix[:,row])
            dissimilarity_matrix[col][row] = 1-correlation(x,y) 
    return dissimilarity_matrix

def shortest_path(edges, f, t):
    """Return tuple of cost and path for the shortest path from node f to t

    edges: list of edge-tuples [(node1,node2,edge_weight)]
    f: string, label of start node
    t: string, label of end node

    This is a python implementation of Dijkstra's shortest path algorithm.
    """
    g = defaultdict(list)
    for l,r,c in edges:
        g[l].append((c,r))

    q, seen = [(0,f,())], set()
    while q:
        (cost,v1,path) = heappop(q)
        if v1 not in seen:
            seen.add(v1)
            path = (v1, path)
            if v1 == t: 
                return (cost, path)
            for c, v2 in g.get(v1, ()):
                if v2 not in seen:
                    heappush(q, (cost+c, v2, path))

    return float("inf")

class Clusterset(object):
    """Clusterset object describing a cluster

    """
    def __init__(self,left=None,right=None,dissimilarity=0.0,ident=None):
        """ 
        ident: identifier of a leaf node (-1 for internal nodes)
        left: clusterset, left child of the node
        right: clusterset, right child of the node
        dissimilarity: dissimilarity between left and right children
        """
        self.left = left
        self.right = right
        self.ident = ident
        self.dissimilarity = dissimilarity

def print_tree(clust,labels=None,n=0):
    """Print graphical representation of clusterset object
    """
    for i in range(n): print ' ',
    if clust.ident<0: # negative id means that this is branch
        print '-'
    else: # positive id means that this is an endpoint
        if labels==None: 
            print clust.ident
        else: 
            print labels[clust.ident]
    # now print the right and left branches
    if clust.left!=None: 
        print_tree(clust.left,labels=labels,n=n+1)
    if clust.right!=None: 
        print_tree(clust.right,labels=labels,n=n+1)

def get_ordered_elements_dissimilarity(clust,list_of_elements=None,\
    list_of_dissim=None):
    """Return ordered list of elements and dissimilarity from clusterset object
    
    clust: Clusterset object
    list_of_elements: list with node identifiers
    list_of_dissim: list of dissimilarity values 
    """
    if list_of_elements is None: 
        list_of_elements=[]
        list_of_dissim=[]
    if clust.ident < 0: # negative id means that this is a branch
        list_of_elements.append('-')
        # append dissimilarity between the right and left brances
        list_of_dissim.append(clust.dissimilarity)
    else: # positive id indicates leaf node (the dissim will be zero)
        list_of_elements.append(clust.ident)
        list_of_dissim.append(clust.dissimilarity)
    if clust.left is not None:
        get_ordered_elements_dissimilarity(clust.left,\
            list_of_elements=list_of_elements,list_of_dissim=list_of_dissim)
    if clust.right is not None: 
        get_ordered_elements_dissimilarity(clust.right,\
            list_of_elements=list_of_elements,list_of_dissim=list_of_dissim)

    return list_of_elements, list_of_dissim

def cut_tree(list_of_elements, list_of_dissim, h=1):
    """Return list of clusters by 'cutting' the dendrogram at height h
    
    list_of_elements: list of node identifiers
    list_of_dissim: list of dissimilarity values
    """
    clustlist = []
    cl = []
    for i in range(len(list_of_elements)):
        if list_of_dissim[i] < h:
            if(list_of_elements[i] != '-'):
                cl.append(list_of_elements[i])
        if list_of_dissim[i] >= h:
            if len(cl) > 0:
                clustlist.append(cl)
            cl=[]
    clustlist.append(cl)  
    return clustlist  

def csv_to_matrix(csv):
    """Creates a matrix from a .csv file""" 

    matrix=[]
    for row in csv:
        matrix.append(row.strip().split(','))
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
           
            matrix[i][j] = float(matrix[i][j])
    
    return matrix
 
def csv_to_edges(lines, weight=1):
    """Return list of edges from csv adjacency matrix
    
    lines: open csv file, which behaves like a list of lines
    
    If node1 and node2 are connected in the adjacency matrix,
    they will be added to the list as (node1,node2,1) and 
    (node2,node1,1). In this way the graph is undirected,
    and all edges have weight 1.
    """
    res = []
    node_i = 0
    for line in lines:
        if not line.strip():
            continue
        parts = map(int,line.strip().split(','))
        for node_j, linked in enumerate(parts):
            if linked:
                res.append((node_i, node_j, weight))
        node_i += 1
    return res
 
def hclust(matrix):
    """Return Clusterset object, representing a hierarchical tree

    matrix: matrix of dissimilarity values from a correlation
    """   
    # initialize the cluster, initially each element is in its own cluster
    clust = [Clusterset(ident=i) for i in range(len(matrix))]
    
    # at each iteration the number of clusters is decreased by one
    # we will only finish when there is a unique element. 

    condition = True
    mindist = 0
    A = None
    B = None

    while len(clust)>1:
        dim = len(matrix)
        
        # step 1: find the two least dissimilar elements, called A and B
        for col in range(dim):
            for row in range(dim):
                if row == 1 and col ==0:
                    mindist = matrix[col][row]
                    A = col
                    B = row
                if matrix.shape[0] > 2:
                    if matrix[col][row] > 0.0:
                        if matrix[col][row] < mindist:
                            mindist = matrix[col][row]
                            A = col
                            B = row
                else:
                    condition = False
                    if matrix[col][row] > 0.0:
                        mindist = matrix[col][row]
                        A = col
                        B = row
                        # step 4 now there is a new cluster formed by elements A and B 
                        # that have dissimilarity mindist and an identifier of -1,
                        # because it is an internal node
                        newcluster = Clusterset(left=clust[A], right=clust[B],\
                        dissimilarity=mindist, ident=-1) 
                        # now delete the old nodes  from the cluster
                        del clust[A]  
                        del clust[B-1]
                        # add the new cluster
                        clust.append(newcluster)
                    break
                    
        if condition:
            #Get complete row with position A
            list_A = list(matrix[A,:])
            #Get complete row with position B
            list_B = list(matrix[B,:])
            dim = len(matrix[A,:])
            new_item = []
            
            # step 2: calculate a new row with the average dissimilarity of 
            # the two joined elements
            for index in range(dim):
                if index != B and index != A:
                    c1 = list_A[index]
                    c2 = list_B[index]
                    c3 = [c1,c2]        

                    new_item.append(round(mean(c3),4))

            # step 3: delete the old rows from the matrix
            matrix = delete(matrix, (B), axis = 0)
            matrix = delete(matrix, (A), axis = 0)
            matrix = delete(matrix, (B), axis = 1)
            matrix = delete(matrix, (A), axis = 1)
            #Inset row at the end
            matrix = vstack((matrix,new_item))
            new_item.append(0.0)
            #Inset column at the end
            matrix = insert(matrix, len(new_item)-1, [new_item], axis = 1)
            # step 4 now there is a new cluster formed by elements A and B 
            # that have dissimilarity mindist and an identifier of -1,
            # because it is an internal node
            newcluster = Clusterset(left=clust[A], right=clust[B],\
            dissimilarity=mindist, ident=-1) 
            # now delete the old nodes  from the cluster
            del clust[A]  
            del clust[B-1]
            # add the new cluster
            clust.append(newcluster)
    return clust[0] 

if __name__ == "__main__": 

    # step 1: go from csv to graph object (already implemented)
    edges = csv_to_edges(open('toy_example.csv'))

    # step 2: go from edges to distance matrix 
    dist_matrix = distance(edges)
    
    # step 3: turn distance matrix to association matrix
    assoc_matrix = association(dist_matrix)
    
    # step 4: calculate dissimilarity matrix from association matrix
    diss_matrix = dissimilarity_matrix_from_correlation(assoc_matrix)
    
    print "##############"
    print "# Question 1 #"
    print "##############"
    print 
    print "Distances"
    print dist_matrix
    print 
    print "Association"
    print assoc_matrix
    print 
    print "Dissimilarities"
    print diss_matrix

    print "##############"
    print "# Question 2 #"
    print "##############"
    print 
    # step 5: hierarchical clustering of the dissimilarity matrix
    clust = hclust(diss_matrix)
    # step 6: visualize/print some clustering results (implemented)
    # possible adjust the parameters to answer the questions
    
    print
    print "Cluster Tree"
    print_tree(clust)
    list_of_elements, list_of_dissim =  get_ordered_elements_dissimilarity(clust) 
    #print list_of_dissim
    print 
    print "height: 0.02 - 0.05"
    print(cut_tree(list_of_elements, list_of_dissim, h=0.02))
    print 
    print "height: 0.14 - 0.18"
    print(cut_tree(list_of_elements, list_of_dissim, h=0.14))
    print 
    print "height: 0.19 - 0.26"
    print(cut_tree(list_of_elements, list_of_dissim, h=0.19))
    print
    print "height: 0.27 >"
    print(cut_tree(list_of_elements, list_of_dissim, h=0.27))

    #print(cut_tree(list_of_elements, list_of_dissim, h=0.15))
    #print(cut_tree(list_of_elements, list_of_dissim, h=0.2))
    
    dist_matrix = csv_to_matrix(open('shortest_path_distance_high_risk_network.csv'))
    csv = open('metabolite_nodes.csv')
    matrix=[]
    labels = {}
    for row in csv:
        matrix.append(row.strip().split(','))
    for i in range(1,len(matrix)):
        for j in range(len(matrix[i])-1):
            idx = int(matrix[i][j].strip('"'))
            val = str(matrix[i][j+1].strip('"'))
            labels[idx] = val

    print "##############"
    print "# Question 3 #"
    print "##############"
    print 
    assoc_matrix = association(dist_matrix)
    diss_matrix = dissimilarity_matrix_from_correlation(assoc_matrix)
    clust = hclust(diss_matrix)
    print_tree(clust,labels)
    list_of_elements, list_of_dissim =  get_ordered_elements_dissimilarity(clust)  
    print
    print "Height: 0.45"
    print(cut_tree(list_of_elements, list_of_dissim, h=0.45))
    print 
    print "Height: 0.3"
    print(cut_tree(list_of_elements, list_of_dissim, h=0.3))
    print
    print "Height: 0.25"
    print(cut_tree(list_of_elements, list_of_dissim, h=0.25))