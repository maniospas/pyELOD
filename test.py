import networkx as nx
import max_conductance
import measures


method = 'ER' #BA or ER
nodes = 100
param = 0.03

elod_fast_conductance = list()
pcst_fast_conductance = list()

print("name = '"+method+"("+str(nodes)+','+str(param)+")'")
for i in range(100):
    if method == 'ER':
        G = nx.erdos_renyi_graph(nodes, param, seed=i).to_directed()
    elif method == 'BA':
        G = nx.generators.barabasi_albert_graph(nodes, param, seed=i).to_directed()
    G.remove_nodes_from(list(nx.isolates(G)))
    G.remove_edges_from(nx.selfloop_edges(G))
    r = list(G.nodes)[0]

    trace, _ = max_conductance.core(G, r, max_conductance.ELODFast())
    elod_fast_conductance.append(measures.conductance(G, trace))
    trace, _ = max_conductance.core(G, r, max_conductance.PCSTFast())
    pcst_fast_conductance.append(measures.conductance(G, trace))

print('ELODFast = ', elod_fast_conductance, ';')
print('PCSTFast = ', pcst_fast_conductance, ';')
print("scatter(PCSTFast, ELODFast, '.')")

print('Errors', sum(1 for i in range(len(elod_fast_conductance)) if elod_fast_conductance[i]<pcst_fast_conductance[i]), '/', len(elod_fast_conductance))
print('Better', sum(1 for i in range(len(elod_fast_conductance)) if elod_fast_conductance[i]>pcst_fast_conductance[i]), '/', len(elod_fast_conductance))

