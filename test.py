import networkx as nx
import max_conductance
import measures
import time


method = 'ER' #BA or ER
nodes = 300
param = 0.03

elod_fast_conductance = list()
pcst_fast_conductance = list()

print("name = '"+method+"("+str(nodes)+','+str(param)+")'")
for i in range(10):
    if method == 'ER':
        G = nx.erdos_renyi_graph(nodes, param, seed=i).to_directed()
    elif method == 'BA':
        G = nx.generators.barabasi_albert_graph(nodes, param, seed=i).to_directed()
    G.remove_nodes_from(list(nx.isolates(G)))
    G.remove_edges_from(nx.selfloop_edges(G))
    r = list(G.nodes)[-1]
    G = max_conductance._remove_cycles(G, r)
    if len(G) < 10:
        continue
    start = time.time()
    trace, _ = max_conductance.core(G, r, max_conductance.ELODFast())
    print(time.time()-start)
    elod_fast_conductance.append(measures.conductance(G, trace))
    start = time.time()
    trace, _ = max_conductance.core(G, r, max_conductance.PCSTFast())
    print(time.time()-start)
    pcst_fast_conductance.append(measures.conductance(G, trace))
    print("----")

elod_fast_conductance = [round(v) for v in elod_fast_conductance]
pcst_fast_conductance = [round(v) for v in pcst_fast_conductance]
print('ELODFast = ', elod_fast_conductance, ';')
print('PCSTFast = ', pcst_fast_conductance, ';')
print("scatter(PCSTFast, ELODFast, '.')")

print('Errors', sum(1 for i in range(len(elod_fast_conductance)) if elod_fast_conductance[i]<pcst_fast_conductance[i]), '/', len(elod_fast_conductance))
print('Better', sum(1 for i in range(len(elod_fast_conductance)) if elod_fast_conductance[i]>pcst_fast_conductance[i]), '/', len(elod_fast_conductance))

