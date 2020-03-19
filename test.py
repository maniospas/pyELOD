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


core_trace_conductance = list()
core_pcst_conductance = list()
core_trace_a = list()
core_pcst_a = list()
core_trace_time = list()
core_pcst_time = list()
trace_conductance = list()
pcst_conductance = list()
trace_time = list()
pcst_time = list()
trace_ELOD = list()
pcst_ELOD = list()
trace2_ELOD = list()
pcst2_ELOD = list()
trace10_ELOD = list()
pcst10_ELOD = list()


for i in range(100):
    #G = nx.erdos_renyi_graph(100, 0.2, seed=i).to_directed()
    G = nx.generators.barabasi_albert_graph(1000, 5, seed=i).to_directed()
    G.remove_nodes_from(list(nx.isolates(G)))
    G.remove_edges_from(nx.selfloop_edges(G))
    r = list(G.nodes)[0]

    tic = time.time()
    T, trace_a = ELOD.core(G, r, ELOD.trace)
    core_trace_time.append(time.time()-tic)
    tic = time.time()
    W, pcst_a = ELOD.core(G, r, ELOD.pcst)
    core_pcst_time.append(time.time()-tic)

    core_trace_conductance.append(ELOD.conductance(G, T))
    core_pcst_conductance.append(ELOD.conductance(G, W))

    core_trace_a.append(trace_a)
    core_pcst_a.append(pcst_a)

    found_a = min(trace_a, pcst_a)

    tic = time.time()
    T = ELOD.trace(G, r, found_a)
    trace_time.append(time.time()-tic)
    tic = time.time()
    W = ELOD.pcst(G, r, found_a)
    pcst_time.append(time.time()-tic)
    

    trace_conductance.append(ELOD.conductance(G, T))
    pcst_conductance.append(ELOD.conductance(G, W))

    trace_ELOD.append(ELOD.ELOD(G, T, found_a))
    pcst_ELOD.append(ELOD.ELOD(G, W, found_a))


    trace2_ELOD.append(ELOD.ELOD(G, ELOD.trace(G, r, 2), 2))
    pcst2_ELOD.append(ELOD.ELOD(G, ELOD.pcst(G, r, 2), 2))

    T = ELOD.trace(G, r, 10)
    W = ELOD.pcst(G, r, 10)
    trace10_ELOD.append(ELOD.ELOD(G, T, 10))
    pcst10_ELOD.append(ELOD.ELOD(G, W, 10))
    #print(set(list(T))-set(list(W)))
    #print(set(list(W))-set(list(T)))

    print(i)

print('Core trace conductance')
print(core_trace_conductance)
print('Core PCST conductance')
print(core_pcst_conductance)

print('Core trace a')
print(core_trace_a)
print('Core PCST a')
print(core_pcst_a)

print('Core trace time')
print(core_trace_time)
print('Core PCST time')
print(core_trace_time)

print('Trace time')
print(trace_time)
print('PCST time')
print(pcst_time)

print('trace conductance')
print(trace_conductance)
print('PCST conductance')
print(pcst_conductance)

print('Trace ELOD')
print(trace_ELOD)
print('PCST ELOD')
print(pcst_ELOD)

print('2Trace ELOD')
print(trace2_ELOD)
print('2PCST ELOD')
print(pcst2_ELOD)

print('10Trace ELOD')
print(trace10_ELOD)
print('10PCST ELOD')
print(pcst10_ELOD)

for i in range(len(trace10_ELOD)):
    if trace10_ELOD[i]<pcst10_ELOD[i]:
        print(i, trace10_ELOD[i], pcst10_ELOD[i])
for i in range(len(trace2_ELOD)):
    if trace2_ELOD[i]<pcst2_ELOD[i]:
        print(2, i, trace2_ELOD[i], pcst2_ELOD[i])
for i in range(len(trace_ELOD)):
    if trace_ELOD[i]<pcst_ELOD[i]:
        print(2, i, trace_ELOD[i], pcst_ELOD[i])







