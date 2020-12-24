def conductance(G, subgraph):
    if subgraph is None:
        return 0
    internal = 0
    outgoing_total = 0
    for v in subgraph:
        internal += subgraph.out_degree[v]
        outgoing_total += G.out_degree[v]
    if internal == 0:
        return 0
    return outgoing_total/internal-1


def _count(G, subgraph):
    internal = 0
    outgoing_total = 0
    for v in subgraph:
        internal += subgraph.out_degree[v]
        outgoing_total += G.out_degree[v]
    return internal, outgoing_total


def ELOD(G, subgraph, a):
    internal = 0
    outgoing_total = 0
    for v in subgraph:
        internal += subgraph.out_degree[v]
        outgoing_total += G.out_degree[v]
    return outgoing_total-a*internal


#evaluate GW https://arxiv.org/pdf/1908.05767.pdf, optimal value 0.7632

