


def prune_leaves(tree, f):
    while tree.count_terminals() > 0:
        altered = False
        for a in tree.get_terminals():
            if f(a):
                tree.prune(a)
                altered = True
        if not altered:
            break


def align_cn_tree(tree, adata):
    prune_leaves(tree, lambda a: a.name not in adata.obs.index)

    cell_ids = []
    for a in tree.get_terminals():
        cell_ids.append(a.name)

    adata = adata[cell_ids]

    return tree, adata
