#!/usr/bin/env python

"""
Algorithms in Bioinformatics
Author: Diego Montiel
Student number: 880505580110
Script to: Genome Assembly
"""
import random

def is_eulerian(graph):
    """
    Function that determines whether a graph is or not Eulerian
    
    Input: Dictionary of a graph nodes and edges
    :return: Boolean True or False 
    """
    count = 0
    condition = True
    for nodes in graph:
        count = len(graph[nodes])
        for edges in graph.values():
            for edge in edges:
                if edge == nodes:
                    count += 1
        if count % 2 != 0:
            condition = False
            return condition
    return condition

def has_eulerian_path(graph):
    """
    Function that determines whether a graph has an Eulerian path
    
    Input: Dictionary of a graph nodes and edges
    :return: Boolean True or False 
    """
    count = 0
    vertice = 0
    condition = True
    for nodes in graph:
        count = len(graph[nodes])
        for edges in graph.values():
            for edge in edges:
                
                if edge == nodes:
                    count += 1
        if count % 2 != 0:
            vertice +=1
        if vertice > 2:
            condition = False
            return condition
    return condition

def get_edges(graph):
    """
    Function that combines key and value of a dictionary
    
    Input: Dictionary of a graph nodes and edges
    :return: List of tuples of the key and value 
    """
    
    paths = []
    for key,value in graph.items():
        for v in value:
            paths.append((key,v))
    return paths

def get_random_value(tuple):
    """
    Function that returns a random choice of a tuple
    
    Input: tuple of integer
    :return: a random tuple choice 
    """
    
    start = random.choice(tuple)
    start = random.choice(start)

    return start

def get_paths(graph):

    """
    Function that give all possible paths of a graph
    
    Input: Dictionary of a graph nodes and edges
    :return: Dictionary of key node and value of edges
    """
    
    graph = get_edges(graph)
    start = get_random_value(graph)
    #start = graph[0][0]
    
    node  = start
    paths = {}
    paths[0] = []
    paths[0].append(start)
    index = 0
    flag  = False    

    while len(graph) > 0:
        for edge in graph:
            if node == edge[0]:
                node = edge[1]
                graph.remove(edge)
                paths[index].append(node)               
                flag = True
                
        if node == start and flag and len(graph) > 0:

                node   = get_random_value(graph)
                #node   = graph[0][0]
                start  = node
                flag   = False
                index += 1
                paths[index] = []
                paths[index].append(node)

    return paths

def find_eulerian_cycle(graph):
  
    """
    Function that returns the cycle path of a graph
    
    Input: Dictionary of a graph nodes and edges
    :return: List of integers of all the path
    """
    
    if is_eulerian(graph):
        paths       = get_paths(graph)
        list_paths  = []
        final_path  = []
        for key, value in paths.items():
            list_paths.append(value)
       
        final_paths = []
        while len(list_paths) > 1:
            for index, value in enumerate(list_paths[0:len(list_paths)-1]):
                for id_path, v in enumerate(value):
                    for i in range(1,len(list_paths)-index):
                        if v in list_paths[index+i]:
                            if list_paths[i].index(v) != 0:
                                index_cut = list_paths[i].index(v)
                                new_list = list_paths[i][index_cut:-1] +\
                                list_paths[i][:index_cut]
                                new_list.append(v)
                                index_list = list_paths[index].index(v)
                                list_paths[index][index_list] = new_list
                                del list_paths[i]
                                break
                                    
                            else:
                                index_list = list_paths[index].index(v)
                                list_paths[index][index_list] = list_paths[i]
                                del list_paths[i]
                                break
        
        list_paths = list_paths[0]
        for path in list_paths:
            if type(path) is list:
                for element in path:
                    final_path.append(element)
            else:
                final_path.append(path)

        return final_path
    else:
        return 'IS NOT EULERIAN, A CYCLE CAN NOT BE FOUND!'

def find_eulerian_path(graph):

    """
    Function that give an eulerian path
    
    Input: Dictionary of a graph nodes and edges
    :return: List of path
    """
    if is_eulerian(graph):
        path = find_eulerian_cycle(graph)
        path = path[1:]
        return path

    elif not has_eulerian_path(graph):
        print 'Graph has not possible solution!'
        return

    return path

def spectrum(st, k):
    """
    Function that returns a graph from spectrum
    
    Input: st = string of characters (DNA)
    Input: k = length of the
    :return: Dictionary of key node and value of edges
    """
    edges = []
    graph = {}
    c = 0
    node = set()
    for i in range(len(st)- k+1):
            edges.append((st[i:i+k-1], st[i+1:i+k]))
            node.add(st[i:i+k-1])
            node.add(st[i+1:i+k])
            for i in node:
                graph[i] = set()
            for e in edges:
                if e[0] in graph:
                    graph[e[0]].add(e[1])
    final_graph = {}
    for k,v in graph.items():
        final_graph[k] = list(v) 
    return final_graph

if __name__ == "__main__":

    # GRAPH FROM FIG 8.22
    #graph_822 = {1:[3],3:[7],10:[9],9:[4],4:[2,5],7:[10,8],2:[1],5:[8],\
    #    8:[4,6],6:[7]}
    graph_822 = {1:[3],3:[7],10:[9],9:[4],4:[2,5],7:[10,8],2:[1],5:[8],\
        8:[4,6],6:[7]}
    # A SLIGHTLY BIGGER GRAPH, NEEDED FOR Q8
    bigger_graph = {'A':['C'],'C':['E'],'F':['G'],'G':['H'],'H':['A','J'],\
        'E':['F','K'],'J':['K'],'K':['H','D'],'D':['E','I'],\
        'I':['B'], 'B':['D']}
    # SPECTRUM FROM FIG 8.20
    s = ['ATG','TGG','TGC','GTG','GGC','GCA','GCG','CGT']
    k = len(s[0])
    s = ''.join(s)
    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file) 
    ##########################
    # Question 1
    ##########################
    print "################"
    print "# Question 1 ###"
    print "################"
    print "Is the Graph of Fig 8.22 Eulerian? "
    print is_eulerian(graph_822)
    print "\n"
    print "################"
    print "# Question 2 ###"
    print "################"
    print "Does the Graph has an Eulerian path? "
    print has_eulerian_path(graph_822)
    print "\n"
    print "################"
    print "# Question 3 ###"
    print "################"
    print "Eulerian path "
    print find_eulerian_path(graph_822)
    print "\n"
    print "################"
    print "# Question 5 ###"
    print "################"
    graph = spectrum(s,k)
    for k, v in graph.items():
        print k, v
    print "\n"
    print "################"
    print "# Question 6 ###"
    print "################"
    print "Is the spectrum graph Eulerian?"
    print is_eulerian(graph)
    print "Has the spectrum graph an Eulerian path?"
    print has_eulerian_path(graph)
    print "\n"
    print "################"
    print "# Question 7 ###"
    print "################"
    print find_eulerian_path(graph)
    print "\n"
    print "################"
    print "# Question 8 ###"
    print "################"
    print "Bigger graph Eulerian Cycle"
    print find_eulerian_cycle(bigger_graph)




