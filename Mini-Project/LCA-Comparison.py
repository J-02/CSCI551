from Bio import Phylo
from Bio import SeqIO
import matplotlib.pyplot as plt
import random
import time
from io import StringIO
import math
from ngesh.random_tree import gen_tree


def timer(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        return result, end_time-start_time

    return wrapper

def create_random_tree(**node):
    """ Create random tree """
    tree = gen_tree(**node, labels='enum')
    newick_str = tree.write(format=1)
    # Convert the Newick string to a Phylo tree
    phylo_tree = Phylo.read(StringIO(newick_str), "newick")
    return phylo_tree
def get_node_names(tree):
    """ Get list of nodes to randomly pair for queries"""
    node_names = [clade.name for clade in tree.find_clades() if clade.is_terminal()]
    return node_names

def generateQueries(tree, queries = 1, rep = 10):
    """Generate tupes of nodes for querying based on the tree, with sample size = rep"""
    names = get_node_names(tree)
    nodeqs = [[select_nodes(names) for i in range(queries)] for r in range(rep)]
    return nodeqs

def select_nodes(nodes, count=2):
    """ Select two nodes randomly for pair for query creation"""
    return random.sample(nodes, count)

def read_fasta(file):
    # reading the fasta seqs
    sequences = {}
    for record in SeqIO.parse(file, "fasta"):
        sequences[record.id] = str(record.seq)
    print(f"{len(sequences)} sequences read:")
    for i in sequences:
        print(i)
    return sequences
def createPhylogenicTree(file):
    """Builds phylogenic tree from .nwk file"""
    tree_path = file
    tree = Phylo.read(tree_path, 'newick')

    # Visualize the tree
    return tree
def assign_names_to_internal_nodes(tree):
    """Assigns names to internal nodes so that LCA's are identifiable"""
    internal_node_count = 1
    for clade in tree.find_clades():
        if not clade.is_terminal() and clade.name is None:
            clade.name = f"Clade-{internal_node_count}"
            internal_node_count += 1

def naiveLCA(tree, queries):
    """Tree traversal based LCA does one pair at a time using paths from root to node it finds the lowest intersection"""
    def find_path(current_node, target_name, path=None):
        """ Finds path from a node to another node"""
        if path is None:
            path = []
        if current_node is None:
            return None
        path.append(current_node)
        if current_node.name == target_name:
            return path
        for next_node in current_node.clades:
            if find_path(next_node, target_name, path):
                return path
        path.pop()
        return None
    def find_lca(root, name1, name2):
        """ Finds intersection of paths"""
        path1 = find_path(root, name1)
        path2 = find_path(root, name2)
        if path1 is None or path2 is None:
            return None
        lca = None
        for node1, node2 in zip(path1, path2):
            if node1 == node2:
                lca = node1
            else:
                break
        return lca
    lcas = [find_lca(tree.root, node1_name, node2_name) for (node1_name, node2_name) in queries]
    return lcas


def tarjansLCA(tree, queries):
    """Uses adjacency lists, Union Find, and depth first search to find and save all lca's """
    # Build adjacency list
    queries = [tuple(query) if isinstance(query, list) else query for query in queries]

    def build_adjacency_list(tree):
        # creating adjaceny list
        adj_list = {}

        def build_list(clade, parent=None):
            # builds list by going from clade to clade and retains parents
            if clade.name not in adj_list:
                adj_list[clade.name] = []
            if parent:
                adj_list[parent.name].append(clade.name)
            for child in clade.clades:
                # recurisive call for children
                build_list(child, clade)

        build_list(tree.root)
        return adj_list

    adj_list = build_adjacency_list(tree)

    # Union-Find structure
    class UnionFind:
        def __init__(self):
            self.parent = {}

        def make_set(self, x):
            self.parent[x] = x

        def find(self, x):
            if self.parent[x] != x:
                self.parent[x] = self.find(self.parent[x])  # Path compression
            return self.parent[x]

        def union(self, x, y):
            self.parent[self.find(x)] = self.find(y)

    uf = UnionFind()
    ancestors = {}
    lca = {}

    # Iterative DFS
    def dfs_iterative(root, adj_list, uf, ancestors, queries):
        stack = [(root, None)]  # Stack for iterative DFS: (node, parent)
        visited = set()

        while stack:
            node, parent = stack.pop()
            if node not in visited:
                visited.add(node)
                uf.make_set(node)
                ancestors[node] = node

                for child in adj_list.get(node, []):
                    stack.append((child, node))

            if parent is not None:
                uf.union(node, parent)
                ancestors[uf.find(node)] = node

        # process the queries with dynamicish program
        for u, v in queries:
            if u in ancestors and v in ancestors:
                lca[(u, v)] = ancestors[uf.find(u)] if uf.find(u) == uf.find(v) else None

    dfs_iterative(tree.root.name, adj_list, uf, ancestors, queries)

    # Return the results of the queries
    return {query: lca[query] for query in queries if query in lca}
def farach_colton_bender_lca(tree, queries):
    # Euler Tour to record the visitation order and levels of nodes
    def euler_tour(node, depth=0):
        if node is None:
            return [], [], {}
        tour, levels, first_occurrences = [node], [depth], {}
        first_occurrences[node.name] = len(tour) - 1
        for child in node.clades:
            t, l, f = euler_tour(child, depth + 1)
            tour.extend(t + [node])
            levels.extend(l + [depth])
            first_occurrences.update(f)
        return tour, levels, first_occurrences

    euler, levels, first_occurrences = euler_tour(tree.root)

    # Build RMQ Structure using Sparse Table
    def build_sparse_table(levels):
        N = len(levels)
        K = math.floor(math.log2(N)) + 1
        st = [[0] * K for _ in range(N)]
        for i in range(N):
            st[i][0] = i
        for j in range(1, K):
            for i in range(N - (1 << j) + 1):
                if levels[st[i][j-1]] < levels[st[i + (1 << (j-1))][j-1]]:
                    st[i][j] = st[i][j-1]
                else:
                    st[i][j] = st[i + (1 << (j-1))][j-1]
        return st

    st = build_sparse_table(levels)

    # LCA Query using RMQ
    def query_lca(u, v):
        l, r = first_occurrences[u], first_occurrences[v]
        if l > r:
            l, r = r, l
        j = int(math.log2(r - l + 1))
        if levels[st[l][j]] < levels[st[r - (1 << j) + 1][j]]:
            return euler[st[l][j]].name
        else:
            return euler[st[r - (1 << j) + 1][j]].name

    # Process queries
    lca_results = {}
    for u, v in queries:
        lca_results[(u, v)] = query_lca(u, v)

    return lca_results


def experiment(name, tree, func, dr):
    tree, queries = tree
    @timer
    def perform(func, pair):
        return func(tree, pair)
    time = 0.0  # Ensure floating-point arithmetic

    for query_list in queries:
        _, t = perform(func, query_list)
        time += t  # Summing total time for all lists

    avgtime = time / sum(len(query_list) for query_list in queries)
    print(f'Function: {func.__name__} \n Avg Time Per Query: {avgtime} \n Nodes: {name} \n Queries: {len(queries[0])} \n DeathRate: {dr}')
    return avgtime

def plot_results(results, variable):
    for algo_name in results:
        x = list(results[algo_name][variable].keys())
        y = list(results[algo_name][variable].values())
        plt.plot(x, y, label=algo_name)

    plt.xlabel(variable.capitalize())
    plt.ylabel('Average Time')
    plt.title(f'Performance of LCA Algorithms with Varying {variable.capitalize()}')
    plt.legend()
    plt.show()

def main():
    # Initialize dictionaries
    treenodes = {}
    treequeries = {}
    results = {}

    # Collect data for varying number of nodes
    start = 10
    stop = 1000
    step = 50
    nodes = start
    while nodes <= stop:
        tree = create_random_tree(num_leaves=nodes, birth=1, death=0.5)
        queries = generateQueries(tree, 10)  # fixed number of queries
        treenodes[nodes] = (tree, queries)
        print(f"generated tree for {nodes} nodes")
        nodes += step - start if nodes == start else step

    # Collect data for varying number of queries
    tree = create_random_tree(num_leaves=1000, birth=1, death=0.5)
    start = 1
    stop = 1000
    step = 50
    num_queries = start
    while num_queries <= stop:
        queries = generateQueries(tree, num_queries)
        treequeries[num_queries] = (tree, queries)
        print(f"generated {num_queries} queries for tree")
        num_queries += step - start if num_queries == start else step

    lcas = [tarjansLCA, farach_colton_bender_lca, naiveLCA]

    # Run experiments and record results
    for lca in lcas:
        algo_name = lca.__name__
        results[algo_name] = {"nodes": {}, "queries": {}, "deathrate": {}}

        # Record results for varying nodes
        for node_count, tree in treenodes.items():
            avg_time = experiment(f"{node_count}", tree, lca, .5)
            results[algo_name]["nodes"][node_count] = avg_time

        # Record results for varying queries
        for query_count, tree in treequeries.items():
            avg_time = experiment("1000", tree, lca, .5)
            results[algo_name]["queries"][query_count] = avg_time

    # Graphing the results
    # Plot for varying nodes
    plot_results(results, "nodes")

    # Plot for varying queries
    plot_results(results, "queries")

    print(results)

