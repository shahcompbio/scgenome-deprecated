from argparse import ArgumentParser
import pandas as pd
from scipy.cluster import hierarchy


def get_args():
    description = 'Generate clone tree on copynumber using hierarchical ' + \
        'agglomerative clustering'
    p = ArgumentParser(description=description)

    p.add_argument('copynumber', help='copynumber tsv file')
    p.add_argument('tree', help='output newick tree file')

    return p.parse_args()


def construct_internal_node_name(node, leaf_names):
    leaf_names = node.pre_order(lambda x: leaf_names[x.id])
    leaf_names.sort()
    return ''.join(leaf_names)


def _tree_to_newick(node, parent_dist, leaf_names, is_root=False):
    if node.is_leaf():
        return '{}:{:.2f}'.format(leaf_names[node.id], parent_dist - node.dist)

    left_newick = _tree_to_newick(node.get_left(), node.dist, leaf_names)
    right_newick = _tree_to_newick(node.get_right(), node.dist, leaf_names)

    node_name = construct_internal_node_name(node, leaf_names)

    if is_root:
        newick = '({},{}){};'.format(left_newick, right_newick, node_name)
    else:
        newick = '({},{}){}:{:.2f}'.format(
            left_newick, right_newick, node_name, parent_dist - node.dist
        )
    return newick


def tree_to_newick(tree, leaf_names):
    return _tree_to_newick(tree, tree.dist, leaf_names, True)


def generate_copynumber_clone_tree(cn):
        # normalize median copynumber
        median_cn = cn.apply(lambda x: x.median(), axis=0)

        X = (cn - median_cn).values.T
        Z = hierarchy.linkage(X, method='complete', metric='cityblock')
        tree = hierarchy.to_tree(Z)
        nw = tree_to_newick(tree, cn.columns)
        return nw


def main():
    argv = get_args()

    cn = pd.read_table(
        argv.copynumber, index_col=['chr', 'start', 'end', 'width'],
        dtype={'chr': str}
    )

    nw = generate_copynumber_clone_tree(cn)

    with open(argv.tree, 'w') as out_f:
        out_f.write(nw)
    return


if __name__ == '__main__':
    main()
