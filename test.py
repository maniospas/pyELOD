import networkx as nx
import GW


def show(G, subgraph=None):
    import matplotlib.pyplot as plt
    nx.draw(G, with_labels=True, node_size=150, width=1.0,
            node_color=[1 if subgraph is None or v in subgraph else 0 for v in G.nodes()],
            edge_color=[1 if subgraph is None or subgraph.has_edge(v,u) else 0 for v,u in G.edges()],
            node_cmap=plt.cm.rainbow,
            edge_cmap=plt.cm.rainbow )
    plt.show()


def MEPT(G, r, a):
    traceELOD = {v: float('-inf') for v in G}
    traceELOD[r] = G.degree[r]
    best_predecessor = {}
    for _,v in nx.bfs_edges(G,r):
        for u in G.predecessors(v):
            if v != u and traceELOD[v] < G.out_degree[v]+traceELOD[u]-a:
                traceELOD[v] = G.out_degree[v]+traceELOD[u]-a
                best_predecessor[v] = u
    T = nx.DiGraph()
    for u, v in best_predecessor.items():
        T.add_edge(v, u)
    return T

def trace(G, r, a):
    T = MEPT(G, r, a)
    traceELOD = {v: G.out_degree[v] for v in G}
    needs_calc = {v: T.out_degree[v] for v in G}
    pending = list([v for v in G if needs_calc[v] == 0])
    while len(pending) != 0:
        v = pending.pop(0)
        preds = list(T.predecessors(v))
        if len(preds)==0:
            continue
        if len(preds) > 1:
            raise Exception("MEPT should be a tree")
        pred = preds[0]
        traceELOD[pred] += max(0, traceELOD[v]-a)
        needs_calc[pred] -= 1
        if needs_calc[pred] == 0:
            pending.append(pred)
    tr = nx.DiGraph()
    pending = [r]
    while len(pending) != 0:
        v = pending.pop(0)
        for u in T.successors(v):
            if traceELOD[u] >= a:
                tr.add_edge(v, u)
                pending.append(u)
    return tr

def least_trace(G, r):
    return trace(G, r, 1.E-12)


def core_trace(G, r):
    eps = 1.E-6
    tr_cond = eps
    tr = None
    while True:
        a = tr_cond + 1
        next_tr = trace(G, r, a)
        if len(next_tr) == 0:
            break
        tr = next_tr
        tr_cond = conductance(G, tr)
    return tr


def GW_trace(G, r):
    tr = None
    best_value = 0
    for _ in range(10):
        next_tr = GW_cut(G, r)
        tr_cond = conductance(G, next_tr)
        if tr_cond > best_value:
            best_value = tr_cond
            tr = next_tr
        print(tr_cond, ' (max', best_value,')')
    return tr


def BFS(G, r):
    T = nx.DiGraph()
    for v,u in nx.bfs_edges(G,r):
        T.add_edge(v, u)
    return T

def ELOD(G, subgraph, a):
    internal = 0
    outgoing_total = 0
    for v in subgraph:
        internal += subgraph.out_degree[v]
        outgoing_total += G.out_degree[v]
    return outgoing_total-a*internal

def GW_cut(G, r):
    return nx.algorithms.traversal.breadth_first_search.bfs_tree(nx.subgraph(G, GW.GW_cut(G, r)), r)

def conductance(G, subgraph):
    internal = 0
    outgoing_total = 0
    for v in subgraph:
        internal += subgraph.out_degree[v]
        outgoing_total += G.out_degree[v]
    if internal == 0:
        return 0
    return outgoing_total/internal-1

#evaluate GW https://arxiv.org/pdf/1908.05767.pdf, optimal value 0.7632


G = nx.les_miserables_graph().to_directed()
r = list(G.nodes)[20]
print(r)
T = core_trace(G, r)
print(ELOD(G,T,2))
print(ELOD(G,BFS(G, r),2))
print('Least Trace Cut Size', ELOD(G, least_trace(G, r), 0))
print('GW Cut Size', ELOD(G, GW_cut(G, r), 0))
print('MEPT Conductance', conductance(G,MEPT(G, r,2)))
print('Trace Conductance', conductance(G,T))
W = GW_cut(G, r)
print('GW Conductance', conductance(G,W))
show(W, None)





