import networkx as nx
from community import community_louvain
from networkx.algorithms import community
import nxmetis
import sys

G = nx.read_edgelist(sys.argv[1], nodetype=int, comments='#')
DG = nx.read_edgelist(sys.argv[1], nodetype=int, comments='#', create_using=nx.DiGraph)

print(nx.is_connected(G))
