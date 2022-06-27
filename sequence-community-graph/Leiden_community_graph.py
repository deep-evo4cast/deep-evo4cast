import argparse
import numpy as np
import igraph as ig
import leidenalg as la

parser =  argparse.ArgumentParser(description='Partition a graph using Leiden algorithm')

parser.add_argument('-n', '--nodes', required=True, help='File containing sequence indices and their ddG values with the L strain sequence (format : Sequence index | ddG value | Node size) ')
parser.add_argument('-e', '--edges', required=True, help='File containing edges with mutational probabilities (format : Sequence 1 index | Sequence 2 index | mutational probability between them)')
parser.add_argument('-lo', '--loi', required=True, help='File containing indices of 59 local optima sequences(PVs) (format : Sequence index)')
parser.add_argument('-pv', '--pvs', required=True, help='File containing indices of 8 active PVs (format : Cluster ID \':\' Sequence index)')
parser.add_argument('-w','--wt', required=True, help='File containing indice of the L strain sequence (format : Cluster ID \':\' Sequence index)' )
parser.add_argument('-seq', '--sequences',required=True, help='File containing all 4507188 PV sequences including the L strain (format : Amino Acid sequence | ddG value)')
parser.add_argument('-t', '--threshold', default=0, required=False, help='Prune edges with defined weight below given threshold')


args = parser.parse_args()

nodes_file = args.nodes
edges_file = args.edges
lo_file = args.loi
pv_file = args.pvs
w_file = args.wt
seq_file = args.sequences
mut_threshold = float(args.threshold)


# %%
##for nodes##
nodes = []
nodes_size =[]


##for edges##
edges = []
edges_proba= []
self_edges= []

# for local optima indices
lo_indices = []
pv_indices = []
w_indices = []

## for sequences
sequences = []


# %%
def local_opt_repartition(lo_list, part):
    lo_repart = []
    for i in range(len(part.sizes())): #number of communities
        lo_repart.append(0)
    for i in lo_list:
        lo_repart[part.membership[i]]+=1
    return lo_repart # number of local optima per community


# %%
g = ig.Graph(directed=True)


with open(nodes_file,'r') as nfile:
    lines=nfile.read().splitlines()
    for line in lines:
        line=line.split()
        nodes.append(line[0])
        nodes_size.append(float(line[2]))
    g.add_vertices(nodes)
    
with open(edges_file,'r') as efile:
    lines=efile.read().splitlines()
    for line in lines:
        line=line.split()
        if line[0]!=line[1] and float(line[2]) > mut_threshold: # delete self edges and prune with threshold
            edges.append((line[0],line[1]))
            edges_proba.append(float(line[2])) # filling edges probabilities list 
    g.add_edges(edges)

if lo_file:
    with open(lo_file, 'r') as lfile:
        lines=lfile.read().splitlines()
        for line in lines:
            lo_indices.append(int(line))

if pv_file:
    with open(pv_file, 'r') as pvfile:
        lines=pvfile.read().splitlines()
        for line in lines:
            pv_indices.append(int(line.split(':')[1]))

if w_file:
    with open(w_file, 'r') as wfile:
        lines=wfile.read().splitlines()
        for line in lines:
            w_indices.append(int(line.split(':')[1]))

            
if seq_file:
    with open(seq_file, 'r') as seqfile:
        lines=seqfile.read().splitlines()
        for line in lines:
            sequences.append(line.split()[0])
g.es["weights"] = edges_proba  # adding weight attributes to the graph


# %%
partition = la.find_partition(g, weights = 'weights', partition_type=la.ModularityVertexPartition, n_iterations=50)

# %%
print('Nb partitions: '+str(len(partition.sizes())))
print(partition.sizes())

# %%
aggregated_partition = partition.aggregate_partition() # one node per community

# %%
part_wt = partition.membership[w_indices[0]]
aggregated_partition.graph.es['color'] = 'light gray'
aggregated_partition.graph.delete_vertices(aggregated_partition.graph.vs.select(_degree_lt=3)) # removing from graph communities with less than 3 neighbors

# %% [markdown]
# Checking how many vertices are left

# %%
print(len(aggregated_partition.graph.vs))

# %% [markdown]
# Community graph, simplified by removing self edges and undirected

# %%
comm_graph = aggregated_partition.graph.simplify().as_undirected()

# %% [markdown]
# Counting edges from WT community to other communities and coloring them in red. Also modifying edge width by adding edge size proportional to log of sum of weights from WT community to others.
# For the communities that are linked to the WT community we check if they contain PVs. If they do, we write a file with all of the sequences in that community.

# %%
nb_out_edges = 0 # number of edges from wt community to other communities
comm_graph.es['edgewidth'] = 0.5
if part_wt < len(comm_graph.vs):
    for i in range(len(comm_graph.vs)):
        if i != part_wt:
            out_edge = comm_graph.es.select(_between=(aggregated_partition[part_wt], aggregated_partition[i]))
            comm_graph.es.select(_between=(aggregated_partition[part_wt], aggregated_partition[i]))['color'] = 'indianred' # coloring edges that connect WT community in red
            comm_graph.es.select(_between=(aggregated_partition[part_wt], aggregated_partition[i]))['edgewidth'] = np.log2(1+aggregated_partition.weight_from_comm(i,part_wt)) # edge size proportional to log of sum of weights from l strain com to others
            if (len(out_edge)==1):
                nb_out_edges+=1
                for j in pv_indices:
                    if part_wt == partition.membership[j]: # looking for the active PV that is in the community of the WT
                        print(j)
                        print(sequences[j])
                    if i == partition.membership[j]: # looking for the active PV that is in the community linked to the WT community
                        print(j)
                        print(sequences[j])
                        with open("pv_community_"+str(i)+".fasta",'w') as f:
                            for k in partition[i]:
                                f.write(">pv"+str(k)+'\n')
                                f.write(sequences[k]+'\n')
            elif (len(out_edge)>1):
                print(len(out_edge))
print("Nb of edges from WT community to other communities : "+str(nb_out_edges))
    

# %%
comm_graph.vs['color']='gainsboro'
comm_graph.vs[part_wt]['color']='red'
comm_graph.vs['framewidth']=0.5
for i in lo_indices:
    if partition.membership[i] < len(comm_graph.vs) and partition.membership[i] != part_wt:
        comm_graph.vs[partition.membership[i]]['color']='lightsteelblue'

for i in pv_indices:
    if partition.membership[i] < len(comm_graph.vs):
        comm_graph.vs[partition.membership[i]]['framewidth']=1.5

# %%
layout = comm_graph.layout_kamada_kawai(maxiter=20000)
ig.plot(comm_graph, "test_w.svg", opacity=0.6, vertex_size = 5*np.log10(partition.sizes()), edge_width=comm_graph.es['edgewidth'], edge_arrow_width=0.5, edge_arrow_size=0.5, edge_color = comm_graph.es['color'], layout = layout, vertex_frame_width=comm_graph.vs['framewidth'])
#ig.plot(comm_graph, opacity=0.6, vertex_size = 5*np.log10(partition.sizes()), edge_width=0.5, edge_arrow_width=0.5, edge_arrow_size=0.5, edge_color = comm_graph.es['color'], layout = layout, vertex_frame_width=comm_graph.vs['framewidth'])



