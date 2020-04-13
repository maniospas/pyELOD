import networkx as nx
import pcst_fast
import numpy as np
from heap import Heap
from measures import conductance


def traversal(G, r):
    pending = [r]
    T = nx.DiGraph()
    T.add_node(r)
    max_dist = len(G)+1
    dist = {v: max_dist for v in G}
    dist[r] = 0
    best_predecessor = {}
    while len(pending) != 0:
        u = pending.pop(0)
        for v in G.successors(u):
            if v not in T:
                pending.append(v)
                T.add_edge(u, v)
                dist[v] = dist[u] + 1
            else:
                preds = list(G.predecessors(v))
                if len(preds) == 0:
                    continue
                if len(preds) > 1:
                    raise Exception("MEPT should be a tree")
                c = preds[0]
                if dist[c] > dist[u]:
                    T.remove_edge(c, v)
                    T.add_edge(u, v)
                    dist[v] = dist[u] + 1

    preds = {v: T.in_degree[v] for v in T}
    pending = [r]
    edges = list()
    while len(pending) != 0:
        u = pending.pop(0)
        for v in T.successors(u):
            preds[v] -= 1
            edges.append((u, v))
            if preds[v] == 0:
                pending.append(v)
    return edges


def MEPT(G, r, a):
    prev_a = a
    a = max(deg for _, deg in G.out_degree())
    subtraceELOD = {v: float('-inf') for v in G}
    subtraceELOD[r] = G.out_degree[r]
    heap = Heap()
    visited = {v: False for v in G}
    subtraceELOD[r] = G.out_degree[r]
    heap.add(r, subtraceELOD[r])
    best_predecessor = {v: list() for v in G}
    for u in heap:
        visited[u] = True
        for v in G.successors(u):
            if not visited[v]:
                subtraceELOD[v] = max(subtraceELOD[v], G.out_degree[v]+subtraceELOD[u]-a)
                heap.add(v, -subtraceELOD[v])
                best_predecessor[v].append(u)

    T = nx.DiGraph()
    for v, u_list in best_predecessor.items():
        for u in u_list:
            T.add_edge(u, v)

    a = prev_a
    subtraceELOD = {v: float('-inf') for v in G}
    subtraceELOD[r] = G.out_degree[r]
    best_predecessor = {}
    pending = [r]
    needs_calc = {v: T.in_degree[v] for v in T}
    needs_calc[r] = 0
    while pending:
        u = pending.pop(0)
        for v in T.successors(u):
            if subtraceELOD[v] < G.out_degree[v]+subtraceELOD[u]-a:
                best_predecessor[v] = u
                subtraceELOD[v] = G.out_degree[v]+subtraceELOD[u]-a
            needs_calc[v] -= 1
            if needs_calc[v] == 0:
                pending.append(v)


    T = nx.DiGraph()
    for v, u in best_predecessor.items():
        T.add_edge(u, v)
    return T


class ELODFastPrevious:
    def __call__(self, G, r, a):
        T = MEPT(G, r, a)
        traceELOD = {v: G.out_degree[v] for v in G}
        needs_calc = {v: T.out_degree[v] if v in T else 0 for v in G}
        pending = list([v for v in G if needs_calc[v] == 0])
        while len(pending) != 0:
            v = pending.pop(0)
            preds = list(T.predecessors(v)) if v in T else list()
            if len(preds) == 0:
                continue
            if len(preds) > 1:
                raise Exception("MEPT should be a tree")
            pred = preds[0]
            traceELOD[pred] += max(0, traceELOD[v]-a)
            needs_calc[pred] -= 1
            if needs_calc[pred] == 0:
                pending.append(pred)
        tr = nx.DiGraph()
        tr.add_node(r)
        pending = [r]
        while len(pending) != 0:
            u = pending.pop(0)
            for v in T.successors(u):
                if traceELOD[v] >= a:
                    tr.add_edge(u, v)
                    pending.append(v)
        return tr


def _backward(T, r, a, D0):
    D = {v: D0.get(v, 0) for v in T}
    need_calc = {v: T.out_degree(v) for v in T}
    pending = [v for v in T if need_calc[v]==0]
    while pending:
        v = pending.pop(0)
        if T.in_degree(v) == 0:
            continue
        if T.in_degree(v) > 1:
            raise Exception("Can only perform the BACKWARD step in trees")
        p = next(T.predecessors(v))
        D[p] += max(D[v]-a, 0)
        need_calc[p] -= 1
        if need_calc[p] == 0:
            pending.append(p)
    return D


def _forward(T, r, a, D):
    trace = nx.DiGraph()
    trace.add_node(r)
    pending = list(T.successors(r))
    while pending:
        v = pending.pop(0)
        if T.in_degree(v) > 1:
            raise Exception("Can only perform the FORWARD step in trees")
        if D[v] >= a:
            trace.add_edge(next(T.predecessors(v)), v)
            pending.extend(T.successors(v))
    return trace


def _trace_tree(T, r, a, D0):
    D = _backward(T, r, a, D0)
    return _forward(T, r, a, D)


def elod_rank(G, r, a):
    D0 = {v: G.out_degree(v) for v in G}
    D = D0
    while True:
        next_D = {v: (D0[v]+sum(max(D[u]-a, 0) for u in G.successors(v)))*0.85+0.15*D0[v] for v in G}
        error = sum(next_D[v]-D[v] for v in G) / len(G)
        D = next_D
        if error < 0.001:
            break
        print(D)
    D = {v: D0[v]+sum(max(D[u]-a, 0) for u in G.successors(v)) for v in G}
    print(D)
    return D


def _super_trace_dijkstra(G, r, a, contribution, prev_super_trace):
    super_trace = nx.DiGraph()
    super_trace.add_node(r)
    pending = Heap()
    pending.add(r, 0)
    D0 = {v: G.out_degree(v) for v in G}
    if contribution is not None:
        D = {v: contribution.get(v, 0) for v in G}
    else:
        D = {v: G.out_degree(v) for v in G}
    for v in pending: # this is heap pop
        for u in G.successors(v):
            if u != r and u not in super_trace:
                super_trace.add_edge(v, u)
                if contribution is not None:
                    D[u] = max(D[u], D0[u] + D[v] - max(contribution[u]-a, 0) - a)
                else:
                    D[u] += D[v] - a
                pending.add(u, -D[u]) # replaces if already there
    return super_trace


class ELODTree:
    def __call__(self, G, r, a):
        D0 = {v: G.out_degree(v) for v in G}
        return _trace_tree(G, r, a, D0)


class ELODFast:
    def __call__(self, G, r, a, max_repeats=10):
        D0 = {v: G.out_degree(v) for v in G}
        D = None
        super_trace = None
        for count in range(max_repeats):
            prev_super_trace = super_trace
            super_trace = _super_trace_dijkstra(G, r, a, D, super_trace)
            D = _backward(super_trace, r, a, D0)
            if prev_super_trace is not None:
                edge_diff = sum(1 for v, u in super_trace.edges() if not prev_super_trace.has_edge(v, u))
                if edge_diff == 0:
                    break
        self.reached_max_repeats = (count == max_repeats-1)
        return _forward(super_trace, r, a, D)


class PCSTFast:
    def __call__(self, G, r, a):
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
        T = nx.Graph()
        for edge in edges_selected:
            T.add_edge(node_map_inv[edges[edge][0]], node_map_inv[edges[edge][1]])
        if len(T) == 0:
            T = nx.DiGraph()
            T.add_node(r)
            return T
        T = nx.traversal.bfs_tree(T, r)
        return T


def convert_to_dag(G, r):
    return nx.traversal.bfs_tree(G, r)


def core(G, r, trace_method, eps=1.E-12):
    tr_cond = -1
    tr = None
    found_a = -1
    while True:
        a = tr_cond+1+eps
        next_tr = trace_method(G, r, a)
        tr_cond = conductance(G, next_tr)
        # print(a, tr_cond)
        if len(next_tr) <= 1:
            break
        tr = next_tr
        found_a = a
    if getattr(trace_method, 'reached_max_repeats', False):
        print('Warning: trace method reached max repetitions in its last run. Results are compromised')

    return tr, found_a