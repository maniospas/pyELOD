import networkx as nx
import time
import ELOD


def show(G, subgraph=None):
    import matplotlib.pyplot as plt
    nx.draw(G, with_labels=True, node_size=150, width=1.0,
            node_color=[1 if subgraph is None or v in subgraph else 0 for v in G.nodes()],
            edge_color=[1 if subgraph is None or subgraph.has_edge(v,u) else 0 for v,u in G.edges()],
            node_cmap=plt.cm.rainbow,
            edge_cmap=plt.cm.rainbow )
    plt.show()


#evaluate GW https://arxiv.org/pdf/1908.05767.pdf, optimal value 0.7632




G = nx.erdos_renyi_graph(100, 0.1, seed=2).to_directed()
r = list(G.nodes)[0]

T, trace_a = ELOD.core(G, r, ELOD.trace)
W, found_a = ELOD.core(G, r, ELOD.pcst)

print('Trace Conductance', ELOD.conductance(G, T))
print('PCST Conductance', ELOD.conductance(G, W))


print('Trace a = ', trace_a)
print('PCST a = ', found_a)

tic = time.time()
T = ELOD.trace(G, r, found_a)
print('Trace time', (time.time()-tic), 's')
tic = time.time()
W = ELOD.pcst(G, r, found_a)
print('PCST Time', (time.time()-tic), 's')

print('Trace ELOD', ELOD.ELOD(G, T, found_a))
print('PCST ELOD', ELOD.ELOD(G, W, found_a))






