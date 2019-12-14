import networkx as nx
import pcst_fast
import numpy as np


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


def core(G, r, trace_method=trace, eps=1.E-6):
    tr_cond = eps
    tr = None
    found_a = None
    while True:
        prev_cond = tr_cond
        a = tr_cond + 1
        next_tr = trace_method(G, r, a)
        if len(next_tr) == 0:
            break
        tr = next_tr
        tr_cond = conductance(G, tr)
        if prev_cond == tr_cond:
            break
        found_a = a
    return tr, found_a


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



def pcst(G, r, a):
    node_map = {}
    node_map_inv = {}
    for v in G:
        node_map_inv[len(node_map)] = v
        node_map[v] = len(node_map)
    max_degree = max(G.out_degree(v) for v in G)
    edges = list()
    for i, j in G.edges():
        edges.append([node_map[i], node_map[j]])
    vertices_selected, edges_selected = pcst_fast.pcst_fast(np.asarray(edges, dtype=np.int64),
                                   np.asarray([G.out_degree(v) for v in G], np.float64),
                                   np.asarray([a for _ in range(len(edges))], np.float64),
                                   node_map[r], 1, 'strong', 0)
    T = nx.DiGraph()
    for edge in edges_selected:
        T.add_edge(node_map_inv[edges[edge][0]], node_map_inv[edges[edge][1]])
    return T#nx.algorithms.traversal.bfs_tree(T, r)


def conductance(G, subgraph):
    internal = 0
    outgoing_total = 0
    for v in subgraph:
        internal += subgraph.out_degree[v]
        outgoing_total += G.out_degree[v]
    if internal == 0:
        return 0
    return outgoing_total/internal-1