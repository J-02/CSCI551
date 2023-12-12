from Bio import SeqIO
from Bio import Phylo
import matplotlib.pyplot as plt

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

    tree_path = file
    tree = Phylo.read(tree_path, 'newick')

    # Visualize the tree
    Phylo.draw(tree)
    return tree


def naiveLCA(tree_file, node1_name, node2_name):
    def find_path(current_node, target_name, path=None):
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


    # Find LCA
    lca = find_lca(tree.root, node1_name, node2_name)
    return lca.name if lca is not None else "No common ancestor found"








read_fasta("sequences.fasta")
tree = createPhylogenicTree("tree.nwk")
node1_name = "Node1_Name"  # Replace with the first node's name
node2_name = "Node2_Name"  # Replace with the second node's name

lca_name = naiveLCA(tree, node1_name, node2_name)

