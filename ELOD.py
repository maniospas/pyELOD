import networkx as nx
import pcst_fast
import numpy as np
import heapq as heapq

class Heap:
    def __init__(self):
        self.pq = []
        self.entry_finder = {}
        self._counter = 0

    def __contains__(self, item):
        return item in self.entry_finder

    def __iter__(self):
        return self

    def __next__(self):
        while self.pq:
            priority, count, item = heapq.heappop(self.pq)
            if item is not None:
                del self.entry_finder[item]
                return item
        raise StopIteration()

    def __len__(self):
        return len(self.pq)

    def remove(self, item):
        entry = self.entry_finder[item]
        entry[-1] = None

    def add(self, item, priority=0):
        if item in self:
            self.remove(item)
        entry = [priority, self._counter, item]
        heapq.heappush(self.pq, entry)
        self.entry_finder[item] = entry
        self._counter += 1



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

    """
    print('ELOD', ELOD(G, T, a))
    subtraceELOD = {v: float('-inf') for v in G}
    subtraceELOD[r] = G.out_degree[r]
    T2 = nx.DiGraph()
    for u, v in nx.traversal.bfs_edges(G, r):
        if subtraceELOD[v] < G.out_degree[v]+subtraceELOD[u]-a:
            subtraceELOD[v] = G.out_degree[v]+subtraceELOD[u]-a
            best_predecessor[v] = u
    for u, v in best_predecessor.items():
        T2.add_edge(v, u)
    print('ELOD', ELOD(G, T2, a))
    """
    return T


def trace(G, r, a):
    T = MEPT(G, r, a)
    traceELOD = {v: G.out_degree[v] for v in G}
    needs_calc = {v: T.out_degree[v] for v in G}
    pending = list([v for v in G if needs_calc[v] == 0])
    while len(pending) != 0:
        v = pending.pop(0)
        preds = list(T.predecessors(v))
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


def core(G, r, trace_method=trace, eps=1.E-6):
    tr_cond = -1
    tr = None
    found_a = None
    while True:
        a = tr_cond+1+eps
        next_tr = trace_method(G, r, a)
        tr_cond = conductance(G, next_tr)
        # print(a, tr_cond)
        if len(next_tr)<=1:
            break
        tr = next_tr
        found_a = a
    return tr, found_a





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
    T = nx.Graph()
    for edge in edges_selected:
        T.add_edge(node_map_inv[edges[edge][0]], node_map_inv[edges[edge][1]])
    if len(T) == 0:
        T = nx.DiGraph()
        T.add_node(r)
        return T
    T = nx.traversal.bfs_tree(T, r)
    return T


def conductance(G, subgraph):
    internal = 0
    outgoing_total = 0
    for v in subgraph:
        internal += subgraph.out_degree[v]
        outgoing_total += G.out_degree[v]
    if internal == 0:
        return 0
    return outgoing_total/internal-1

def count(G, subgraph):
    internal = 0
    outgoing_total = 0
    for v in subgraph:
        internal += subgraph.out_degree[v]
        outgoing_total += G.out_degree[v]
    return internal, outgoing_total