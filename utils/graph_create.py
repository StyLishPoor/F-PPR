import networkx as nx
from community import community_louvain
from networkx.algorithms import community
import nxmetis
import sys

G = nx.read_edgelist(sys.argv[1], nodetype=int, comments='#')
DG = nx.read_edgelist(sys.argv[1], nodetype=int, comments='#', create_using=nx.DiGraph)
#DG = G.to_directed()
all_edges = set(DG.edges)
frontier_edges = set(DG.edges) # frontiers edges
louvain_list = []

# partition by louvain
# louvain_dict = community_louvain.best_partition(G, resolution=2)
# max_comm = max(louvain_dict.values()) + 1
# for comm in range(max_comm):
#   louvain_list.append([k for k, v in louvain_dict.items() if v == comm])
#   print(len(louvain_list[comm]))


# partition by metis
_, louvain_list = nxmetis.partition(G, int(sys.argv[2]))

# partition by girvan-newman
# comp = community.girvan_newman(G)
# louvain_list = tuple(sorted(c) for c in next(comp))
# print(len(louvain_list))

partition_num = len(louvain_list)
for comm in range(partition_num):
  print(len(louvain_list[comm]))

GM_id = {} # GM_id[v]: GM's id of vertex v
louvain_internal_edges = [] # louvain_internal_edges[i]: internal edges of GM i
for i, part in enumerate(louvain_list):
  subG = DG.subgraph(part)
  louvain_internal_edges.append(list(subG.edges))
  frontier_edges -= set(subG.edges)
  for v in part:
    GM_id[v] = i

frontier_id = {} # frontier_id[x][y]: edges from GM x to GM y
for i in range(partition_num):
  frontier_id[i] = {}
  for j in range(partition_num):
    if j != i:
      frontier_id[i][j] = []

frontier_nodes = set()
for (inv, outv) in frontier_edges:
  frontier_id[GM_id[inv]][GM_id[outv]].append((inv,outv))
  frontier_nodes.add(inv)
  frontier_nodes.add(outv)
  
for i in range(partition_num):
  f = open(sys.argv[3] + "-" + str(i), 'w')
  for (inv, outv) in louvain_internal_edges[i]:
    f.write(str(inv) + " " + str(outv) + '\n')
  f.write("-1 -1\n")
  for GM_id, edges in frontier_id[i].items():
    for (inv, outv) in edges:
      f.write(str(GM_id) + " " + str(inv) + " " + str(outv) + '\n')
  f.write("-1 -1 -1\n")
  gm_frontier = set(louvain_list[i]) & frontier_nodes
  for fr in gm_frontier:
    f.write(str(fr)+'\n')

