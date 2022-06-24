import igraph as ig
import leidenalg as la
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import random
import numpy as np
from sklearn.manifold import TSNE
import argparse

parser = argparse.ArgumentParser(description = 'Plot T-SNE cluster map')
parser.add_argument('-i', '--input', required=True, help = 'Cluster of each local minimum')
parser.add_argument('-l', '--lon', required=True, help = 'Input file containing all local optima RBD interface residues (format : sequences | energy)')
parser.add_argument('-c', '--clusters', required=True, help = 'Input file containing all cluster representatives (format : line ID (number-1) in lon | seq | ddg)')
parser.add_argument('-g', '--goodclust', required=True, help = 'Input file containing IDs of the cluster representatives that worked (format : Cluster ID)')
parser.add_argument('-w', '--wildtype', required=True, help = 'Input file containing L strain RBD interface residues (format : sequence | energy)')

args = parser.parse_args()

lonfile = args.lon
cfile = args.clusters
gfile = args.goodclust
wfile = args.wildtype
afile = args.input

def hamming(a,b):
    cpt = 0
    for i in range(0,len(a)-1):
        if (b[i]!=a[i]):
            cpt+=1
    if (a[len(a)-1]!=0) and (a[len(a)-1]!=0):
        if (a[len(a)-1]!=b[len(b)-1]):
            cpt+=20
    return cpt

def true_hamming(a,b):
    cpt = 0
    for i in range(0,len(a)-1):
        if (b[i]!=a[i]):
            cpt+=1
    return cpt


tsne = TSNE(perplexity = 6, early_exaggeration = 12, metric=hamming, square_distances = True, learning_rate = 50, n_iter = 1000)

AAs = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8, 'I':9, 'L':10, 'K':11, 'M':12, 'F':13, 'P':14, 'S':15, 'T':16, 'W':17, 'Y':18, 'V':19}


    
sequences = []
with open(lonfile,'r') as f:
    lines = f.read().splitlines()
    for line in lines:
        sequences.append(line.split()[0])

seq2clust = {}
clustIDs = []
clust_sizes = {}
with open(afile,'r') as f:
    lines = f.read().splitlines()
    for line in lines:
        line = line.split()
        seq2clust[int(line[0])] = int(line[1])
        clustIDs.append(int(line[1]))
        if line[1] in clust_sizes.keys():
            clust_sizes[line[1]]+=1
        else:
            clust_sizes[line[1]]=1
            
    clustIDs.append(0)

cIDs = []
with open(cfile,'r') as f:
    lines = f.read().splitlines()
    for line in lines:
        cIDs.append(int(line.split()[0]))

gIDs = []
with open(gfile,'r') as f:
    lines = f.read().splitlines()
    for line in lines:
        gIDs.append(int(line))

    
sorted_clusts = dict(sorted(clust_sizes.items(), key=lambda item: item[1], reverse = True))


cluster_names = {}
cpt = 1
for item in sorted_clusts.items():
    cluster_names[item[0]] = 'PV'+str(cpt)
    cpt+=1


cluster_labels = []
for i in range(0, len(cluster_names)):
    cluster_labels.append(cluster_names[str(clustIDs[cIDs[i]])])

for i,l in enumerate(cluster_labels):
        print(str(i+1)+'->'+str(l))
wt_sequence = []
with open(wfile,'r') as f:
    lines = f.read().splitlines()
    for line in lines:
        wt_sequence = line.split()[0]

seq2vec = []
for i,s in enumerate(sequences):
    v = []
    for a in s:
        v.append(AAs[a])
    v.append(seq2clust[i])
    seq2vec.append(v)

v = []
for a in wt_sequence:
    v.append(AAs[a])
v.append(0)
seq2vec.append(v)

    
ts_embedding = tsne.fit_transform(seq2vec)

clusters_x = []
clusters_y = []
for c in cIDs:
    clusters_x.append(np.array(ts_embedding[c,0]))
    clusters_y.append(np.array(ts_embedding[c,1]))

good_x = []
good_y = []
for g in gIDs:
    good_x.append(np.array(ts_embedding[cIDs[g],0]))
    good_y.append(np.array(ts_embedding[cIDs[g],1]))

nb_seq = len(ts_embedding)

fig, ax = plt.subplots()
mycycler = plt.cycler("color", plt.cm.tab20b.colors)

colors = []
for i in mycycler:
    colors.append(i['color'])
mycolors = []
for i in range(0,59):
    mycolors.append(colors[i%20])
    
mycmap = LinearSegmentedColormap.from_list('my list', mycolors, N = 59)

sc = plt.scatter(ts_embedding[:,0],ts_embedding[:,1], c = clustIDs, cmap = mycmap, alpha = 0.3, s=20, linewidth=0)
ax.scatter(ts_embedding[nb_seq-1,0],ts_embedding[nb_seq-1,1], color = 'black', label = 'L strain variant (L)', alpha = 1, s=30)
ax.scatter(clusters_x,clusters_y, color = 'firebrick', alpha = 1, s=30, label = 'Potential variants (PV)')
for i, txt in enumerate(cluster_labels):
    ax.annotate(txt,(clusters_x[i], clusters_y[i]), size=9)
    
ax.scatter(good_x,good_y, color = 'limegreen', alpha = 1, s=30, label = 'Potential and active variants (PV)')
ax.annotate('L', xy = (ts_embedding[nb_seq-1,0],ts_embedding[nb_seq-1,1]), xytext=(0,-10), textcoords = 'offset pixels')

plt.legend(loc='upper left', fontsize = 'small')
plt.axis('off')
plt.savefig('clustermap.svg')




