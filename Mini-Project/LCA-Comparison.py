from Bio import Phylo
from Bio import SeqIO
import matplotlib.pyplot as plt
import random
import time
from io import StringIO
import math
from ngesh.random_tree import gen_tree
from tralda.datastructures import Tree
from tralda.datastructures import LCA



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
    return {(q[0], q[1]): lca.name for q, lca in zip(queries, lcas) if lca is not None}

class UnionFind:
    """A simple Union-Find class for Tarjan's LCA algorithm."""
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def make_set(self, x):
        self.parent[x] = x
        self.rank[x] = 0

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        root_x = self.find(x)
        root_y = self.find(y)

        if root_x != root_y:
            if self.rank[root_x] > self.rank[root_y]:
                self.parent[root_y] = root_x
            elif self.rank[root_x] < self.rank[root_y]:
                self.parent[root_x] = root_y
            else:
                self.parent[root_y] = root_x
                self.rank[root_x] += 1

def tarjansLCA(tree, queries):
    ancestor = {}
    union_find = UnionFind()
    lca_results = {}
    visited = set()

    def dfs(node):
        # Adjust node identification method here
        node_id = node.name  # or another unique identifier

        union_find.make_set(node_id)
        ancestor[union_find.find(node_id)] = node_id

        # Adjust child node traversal based on Phylo tree structure
        for child in node.clades:
            dfs(child)
            union_find.union(node_id, child.name)
            ancestor[union_find.find(node_id)] = node_id

        visited.add(node_id)

        for (node1_name, node2_name) in queries:
            if node_id in (node1_name, node2_name):
                other = node1_name if node_id == node2_name else node2_name
                if other in visited:
                    lca_results[(node1_name, node2_name)] = ancestor[union_find.find(other)]

    # Start DFS from the root of the Phylo tree
    dfs(tree.root)
    return lca_results

def farachLCA(tree, queries):
    lca = LCA(tree)
    lcas = {(n1, n2): lca.get(n1,n2).label for (n1, n2) in queries}
    return lcas

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
    step = 250
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
    step = 250
    num_queries = start
    while num_queries <= stop:
        queries = generateQueries(tree, num_queries)
        treequeries[num_queries] = (tree, queries)
        print(f"generated {num_queries} queries for tree")
        num_queries += step - start if num_queries == start else step

    lcas = [farachLCA, tarjansLCA, naiveLCA]

    # Run experiments and record results
    for lca in lcas:
        algo_name = lca.__name__
        results[algo_name] = {"nodes": {}, "queries": {}, "deathrate": {}}

        # Record results for varying nodes
        for node_count, tree in treenodes.items():
            if lca == farachLCA:
                handle = StringIO()
                Phylo.write(tree, handle, "newick")
                newick_string = handle.getvalue()
                tree = Tree.parse_newick(newick_string)
            avg_time = experiment(f"{node_count}", tree, lca, .5)
            results[algo_name]["nodes"][node_count] = avg_time

        # Record results for varying queries
        for query_count, tree in treequeries.items():
            if lca == farachLCA:
                handle = StringIO()
                Phylo.write(tree, handle, "newick")
                newick_string = handle.getvalue()
                tree = Tree.parse_newick(newick_string)
            avg_time = experiment("1000", tree, lca, .5)
            results[algo_name]["queries"][query_count] = avg_time

    # Graphing the results
    # Plot for varying nodes
    plot_results(results, "nodes")

    # Plot for varying queries
    plot_results(results, "queries")

    print(results)

main()

# testing all methods return the same answers
tree = createPhylogenicTree("tree.nwk")
q = generateQueries(tree, 10)
assign_names_to_internal_nodes(tree)
for clade in tree.find_clades():
    clade.confidence = None
handle = StringIO()
Phylo.write(tree, handle, "newick")
newick_string = handle.getvalue()
tree1 = Tree.parse_newick(newick_string)

lca_name = naiveLCA(tree, q[1])
lca_name1 = farachLCA(tree1, q[1])
lca_name2 = tarjansLCA(tree, q[1])

for (i,ii) in q[1]:
    print(f"LCA for {i, ii}: {lca_name[(i,ii)]} ")

Phylo.draw(tree)
Phylo.draw(tree1)

