import fminmax as m
import unittest
import numpy as np
import higra as hg
import time

def alpha(cur_marker):
    func = (1+ 1e-9)/(cur_marker + 1e-9)
    return func

def fminmax_indx(tree, altitudes, marker):
    """
    Calculate the fminmax distance between a tree and a marker, and returns
    the fminmax distances and the closest vertex for all nodes on the tree.
    Inputs: 'tree': a binary partition tree
            'altitudes': vector with altitudes for the tree nodes
            'marker': a vector with values from 0 to 1 for each leaf of the tree
    Outputs: 'output': fminmax distances
             'output_idx': closest vertex x* for each node
    """
    output = np.zeros(tree.num_vertices())
    output[tree.root()] = np.inf
    output_idx = np.zeros(tree.num_vertices())
    #idx of min on leaves is the leaf itself
    min_alpha_idx = np.ones(tree.num_vertices())*-1
    min_alpha_idx[:tree.num_leaves()] = range(tree.num_leaves())
    # transform marker and comput min in each region (min(alpha(f(x)))
    alpha_marker = alpha(marker)
    min_alpha_marker = hg.accumulate_sequential(tree, alpha_marker, hg.Accumulators.min)
    max_mu = hg.accumulate_sequential(tree, marker, hg.Accumulators.max)
    # sibling of each node
    sibling = hg.attribute_sibling(tree)
    #parent of each node
    parents = tree.parents()

    #get the leaf that minimizes alpha(f(x)) on each region
    for i in tree.leaves_to_root_iterator(include_leaves = False):
        if(min_alpha_marker[tree.child(0, i)] <= min_alpha_marker[tree.child(1, i)]):
            min_alpha_idx[i] = min_alpha_idx[tree.child(0, i)]
        else:
            min_alpha_idx[i] = min_alpha_idx[tree.child(1, i)]
    output_idx[tree.root()] = min_alpha_idx[tree.root()]
    
    #Calculating distance from node to every other node
    for i in tree.root_to_leaves_iterator(include_root=False):
        # distance to sibling fminmax( marker | sib(n), n) = alpha(maxf(sibling(n))) * (1 + altitudes(parent(n)))
        dist_sib = min_alpha_marker[sibling[i]]*(1 - max_mu[sibling[i]] + altitudes[parents[i]])
        if (output[parents[i]] <= dist_sib):
            output_idx[i] = output_idx[parents[i]]
            output[i] = output[parents[i]]
        else:
            output_idx[i] = min_alpha_idx[sibling[i]]
            output[i] = dist_sib
    #Check if the distance from leaves to the rest of the tree
    #are not bigger than the distance from the leaf to itself
    for i in range(tree.num_leaves()):
        dist_to_itself = min_alpha_marker[i]*(1-marker[i])
        if(dist_to_itself <= output[i]):
            output[i] = dist_to_itself
            output_idx[i] = i

    return output,output_idx.astype(int)


for i in range(10):
    s = 200
    num_seeds = 60
    g = hg.get_4_adjacency_graph((s,s))
    tree, altitudes = hg.bpt_canonical(g, np.random.rand(g.num_edges()))
    markers = np.random.rand(g.num_vertices())
    non_zero_markers = np.random.choice(np.arange(g.num_vertices()), num_seeds, replace=False)
    non_zero_obj = non_zero_markers[:num_seeds//2]
    mo = np.zeros_like(markers)
    mo[non_zero_obj] = markers[non_zero_obj]
    begin_time = time.time()
    py_version = fminmax_indx(tree, np.array(altitudes), mo)[1]
    mid_time = time.time()
    alpha_mo = alpha(mo)
    cpp_return = m.fminmax_idx(tree, np.array(altitudes), np.array(mo), np.array(alpha_mo))

    closest_vertex = cpp_return[0]
    pass_edge_w = cpp_return[1]
    end_time = time.time()
    
    pass_edges = np.zeros(tree.num_vertices())
    lca = hg.make_lca_fast(tree)
    itr = list(range(len(pass_edges)))
    lcas = lca.lca(itr,closest_vertex)
    pass_edges = np.array(lcas).astype(int)
    print(np.allclose(pass_edge_w,altitudes[pass_edges]))
    print(np.allclose(py_version,closest_vertex))
    print("Time for py version: ",str(mid_time-begin_time))
    print("Time for cpp version: ",str(end_time-mid_time))
    print("Speedup = ", (mid_time-begin_time)/(end_time-mid_time))