import networkx as nx
import sys

G = nx.read_edgelist(sys.argv[1], nodetype=int, create_using=nx.DiGraph)

for (u, v) in list(G.edges):
  print(u, v)
