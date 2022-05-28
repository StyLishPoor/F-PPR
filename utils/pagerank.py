import networkx as nx
import sys
G = nx.read_edgelist(sys.argv[1], nodetype=int, comments='#', create_using=nx.DiGraph)
#G = nx.read_edgelist(sys.argv[1], nodetype=int, comments='#')
source = int(sys.argv[2])
pr = nx.pagerank(G, alpha=0.8, nstart={source:1}, personalization={source:1}, tol=1e-08,  dangling={source:1})

for v, val in pr.items():
  print(v, val)
