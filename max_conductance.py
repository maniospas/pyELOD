import networkx as nx
import pcst_fast
import numpy as np
from heap import Heap
from measures import conductance
from collections import deque as Que


def _backward(T, a, D0):
    D = {v: D0.get(v, 0) for v in T}
    need_calc = {v: T.out_degree(v) for v in T}
    pending = [v for v in T if need_calc[v] == 0]
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


def _remove_cycles(G, r):
    a = max(deg for _, deg in G.out_degree())
    subtraceELOD = {v: float('-inf') for v in G}
    subtraceELOD[r] = G.out_degree[r]
    heap = Heap()
    visited = {v: False for v in G}
    subtraceELOD[r] = G.out_degree[r]
    heap.add(r, subtraceELOD[r])
    T = nx.DiGraph()
    T.add_node(r)
    for u in heap:
        visited[u] = True
        for v in G.successors(u):
            if not visited[v]:
                subtraceELOD[v] = max(subtraceELOD[v], G.out_degree[v]+subtraceELOD[u]-a)
                heap.add(v, -subtraceELOD[v])
                T.add_edge(u, v)
    return T


def _traverse_order(G, r, D0):
    in_degree = {v: G.in_degree[v] for v in G}
    should_visit = {v: in_degree[v] for v in G}
    pending = Que([r])
    pending_multi_in = Que()
    visit_order = list()
    while True:
        if pending:
            u = pending.popleft()
        elif pending_multi_in:
            u = pending_multi_in.popleft()
        else:
            break
        visit_order.append(u)
        for v in G.successors(u):
            should_visit[v] -= 1
            if should_visit[v] == 0:
                if in_degree[v] > 1:
                    pending_multi_in.append(v)
                else:
                    pending.append(v)
    if len(visit_order) != len(G):
        raise Exception("Topological order is defined only for connected graphs")
    return visit_order


class ELODFast:
    def __call__(self, G, r, a):
        D0 = {v: G.out_degree[v] for v in G}
        G = _remove_cycles(G, r)
        self.D0 = D0
        self.G = G
        self.a = a
        return self._find_trace(r)

    def _find_trace(self, r):
        T = nx.DiGraph()
        pending = [r]
        for u in pending:
            for v in self.G.successors(u):
                T.add_edge(u, v)
                pending.append(v)


def _trace_tree(T, r, a, D0):
    D = _backward(T, a, D0)
    return _forward(T, r, a, D)


class ELODTree:
    def __call__(self, G, r, a):
        D0 = {v: G.out_degree(v) for v in G}
        return _trace_tree(G, r, a, D0)

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