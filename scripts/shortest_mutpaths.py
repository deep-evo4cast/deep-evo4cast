import argparse
from math import exp
from heapq import *
import numpy as np
from scipy.linalg import expm

parser = argparse.ArgumentParser(description = 'Prepare for Dijkstra')
parser.add_argument('-i', '--input', required=True, help = 'Input file with energies and neighbors IDs (format : ddG value | neighbors)')
parser.add_argument('-e', '--enumeration', required=True, help = 'Input file with sequences enumerated (format : Amino Acid sequence | ddG value)')
parser.add_argument('-m', '--mutprobas', required=True, help = 'File containing mutational probabilities (converted into energies) (format : Mutational probability converted into energy of a sequence to each of its neighbors)')
parser.add_argument('-s', '--start', required=True, help = 'Starting point ID (format : Cluster ID \':\' Sequence index)')
parser.add_argument('-g', '--goal', required=True, help = 'Goal ID (format : Cluster ID \':\' Sequence index)')
parser.add_argument('-r', '--results', required=True, help = 'Output file')

#parser.add_argument('-e', '--energies', required=True, help = '')
args = parser.parse_args()

neighbors_file = args.input
start = args.start
goal = args.goal
mutprobasfile = args.mutprobas
results = args.results
energies = []
sequences = []
AA_order = "IMTNKSRLPHQVADEGFYCW_" # codons order with STOP at the end

with open(args.enumeration, 'r') as f:
    lines = f.read().splitlines()
    for line in lines:
        sequences.append(line.split()[0])
        

a2a_energy = np.loadtxt(mutprobasfile)

# Find the mutation proba "energy" for two neighbor sequences
def getA2Aenergy(s1, s2):
    seq1 = sequences[s1]
    seq2 = sequences[s2]
    findmut = False
    for i1, aa1 in enumerate(seq1):
        if aa1!=seq2[i1]:
            if not findmut:
                mut1 = aa1
                mut2 = seq2[i1]
                findmut = True
            else:
                print("WARNING: there is more than one mutation between neighbors !")
    for i, aa in enumerate(AA_order):
        if aa == mut1:
            imut1 = i
        if aa == mut2:
            imut2 = i
    return a2a_energy[imut1, imut2]


def create_graph(neighbors_file):
    cpt = 0
    neighbors = []
    
    with open(neighbors_file,'r') as nfile:
        lines = nfile.read().splitlines()
        for line in lines:
            line = line.split()
            intline = list(map(int,line[1:]))
            neighbors.append(intline) # for each sequence, store distance/ID tuple of each neighbor in a dictionnary
            energies.append(float(line[0]))
            cpt+=1


        dis_neigh = {}

        for cpt,neigh in enumerate(neighbors):
            dis_neigh[cpt] = [(getA2Aenergy(cpt,x), x) for x in neigh]
        return(dis_neigh)

# Shortest path using Dijkstra algorithm

def dijkstra (u, v):
    visited = set()
    d = {u: 0}
    p = {}
    candidates = [(0, u)] 

    while candidates != []:

        dx, x = heappop(candidates) # retrieve the root of the heap
        if x in visited:
            continue

        visited.add(x)

        for w, y in neighbors(x):
            if y in visited:
                continue
            dy = dx + w
            if y not in d or d[y] > dy:
                d[y] = dy
                heappush(candidates, (dy, y)) # place (or update) y in the heap
                p[y] = x

    path = [v]
    x = v
    while x != u:
        x = p[x]
        path.insert(0, x)
    
    return d[v], path



def neighbors(v):
    return graph[v]


graph = create_graph(neighbors_file)


startpoints = {}
with open(start,'r') as sf:
    lines = sf.read().splitlines()
    for l in lines:
        l = l.split(":")
        startpoints[l[0]] = l[1]

endpoints = {}
with open(goal,'r') as sf:
    lines = sf.read().splitlines()
    for l in lines:
        l = l.split(":")
        endpoints[l[0]] = l[1]


with open(results,'w') as f:
    for s in startpoints:
        for g in endpoints:
            print(s)
            print(g)
            distance, path = dijkstra(int(startpoints[s]),int(endpoints[g]))
            print(distance)
            best_dis_neigh = [(energies[x],x) for x in path]
            f.write("start: "+s+" end: cluster "+g+" distance: "+str(distance)+"\n")
            for (e,i) in best_dis_neigh:
                f.write(str(e)+" "+str(i)+"\n")
   
