


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
    pending = Heap()
    pending.add(r, 0)
    D0 = {v: G.out_degree(v) for v in G}
    D = {v: D0[v] for v in G}
    D[r] = D0[r]
    coming_from = dict()
    for v in pending: # this is heap pop
        if v not in super_trace:
            super_trace.add_node(v)
            if v != r:
                super_trace.add_edge(coming_from[v], v)
            for u in G.successors(v):
                if u not in super_trace:
                    prev_D = D[u]
                    if contribution is not None:
                        D[u] = max(D[u], D0[u] + D[v] + contribution[v]-D0[v] - a - max(contribution[u]-a, 0))
                    else:
                        D[u] = max(D[u], D0[u] + D[v] - a)
                    if D[u] != prev_D or u not in coming_from:
                        coming_from[u] = v
                    pending.add(u, -D[u]) # replaces if already there
    return super_trace



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




def _backward_dag(G, a, D0):
    D = {v: D0.get(v, 0) for v in G}
    need_calc = {v: G.out_degree(v) for v in G}
    pending = [v for v in G if need_calc[v] == 0]
    accounted_successors = {v: dict() for v in G}
    #accounted_times = {v: G.in_degree[v] for v in G}
    while pending:
        v = pending.pop(0)
        for p in G.predecessors(v):
            Dv = D[v]
            for s in accounted_successors[v]:
                if s in accounted_successors[p]:
                    Dv -= accounted_successors[v][s] - a
            D[v] = Dv
            if Dv >= a:
                for s in accounted_successors[v]:
                    if s not in accounted_successors[p]:
                        accounted_successors[p][s] = accounted_successors[v][s]
                        if Dv - accounted_successors[v][s] < a:
                            accounted_successors[p][s] = Dv
                if G.in_degree[v] > 1:
                    accounted_successors[p].add(v)
                D[p] += Dv - a
            need_calc[p] -= 1
            if need_calc[p] == 0:
                pending.append(p)
    return D


def _forward_dag(G, r, a, D):
    trace = nx.DiGraph()
    trace.add_node(r)
    pending = Heap()
    pending.add(r, 0)
    visited = {v: False for v in D}
    for u in pending:
        for v in G.successors(u):
            if D[v] >= a and not visited[v]:
                pending.add(v, -D[v])
                trace.add_edge(u, v)
                visited[v] = True
    print(D[r], ELOD(G, trace, a))
    return trace



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
        print(count)
        return _forward(super_trace, r, a, D)