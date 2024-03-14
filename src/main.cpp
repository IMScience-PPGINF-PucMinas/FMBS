#include "pybind11/pybind11.h"
#include "higra/graph.hpp"
#include "higra/accumulator/accumulator.hpp"
#include "higra/accumulator/tree_accumulator.hpp"
#include "higra/attribute/tree_attribute.hpp"
#include <limits>

#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"

namespace py = pybind11;


auto fminmax_idx(const hg::tree & tree, const xt::pyarray<double> & altitudes, const xt::pyarray<double> & marker,const xt::pyarray<double> & alpha_marker){
    /*
    Calculate the fminmax distance between a tree and a marker, and returns
    the fminmax distances and the closest vertex for all nodes on the tree.
    Inputs: 'tree': a binary partition tree
            'altitudes': vector with altitudes for the tree nodes
            'marker': a vector with values from 0 to 1 for each leaf of the tree
            'alpha_marker': vector with values for the function alpha(marker)
    Outputs:'closest_vertex': closest vertex x* for each node
            'pass_edge_weight': weight of the pass edge between node i and closest vertex i*
    
    */

	xt::pyarray<double> fmm_distances = xt::pyarray<double>::from_shape({hg::num_vertices(tree)});
	xt::pyarray<hg::index_t> closest_vertex = xt::pyarray<hg::index_t>::from_shape({hg::num_vertices(tree)});

    xt::pyarray<double> pass_edge_weight = xt::pyarray<double>::from_shape({hg::num_vertices(tree)});
	fmm_distances[hg::root(tree)] = std::numeric_limits<double>::infinity();

    //When using the code below everything works fine.
    xt::pyarray<hg::index_t> min_alpha_idx = xt::ones<hg::index_t>({hg::num_vertices(tree)});
    min_alpha_idx *= -1;
    for(hg::index_t i = 0;i<hg::num_leaves(tree);i++){
        min_alpha_idx[i] = i;
    }
	//xt::pyarray<hg::index_t> min_alpha_idx = xt::pyarray<hg::index_t>::from_shape({hg::num_vertices(tree)});
    //The following 2 lined return an error on my code
	//xt::view(min_alpha_idx, xt::range(hg::num_leaves(tree), hg::num_vertices(tree)) = -1;
    //xt::view(min_alpha_idx, xt::range(0, hg::num_leaves(tree)) = xt::arange<hg::index_t>(0, hg::num_leaves(tree));

    //transform marker and comput min in each region (min(alpha(f(x)))
    //xt::pyarray<double> alpha_marker = (1 + 1e-9)/(marker + 1e-9);
    auto max_mu = hg::accumulate_sequential(tree, marker, hg::accumulator_max());
    auto min_alpha_marker = hg::accumulate_sequential(tree, alpha_marker, hg::accumulator_min());
    auto sibling = hg::attribute_sibling(tree);
    //Compute children explicitly
    //tree.compute_children();
	auto & parents = tree.parents();
    //calculate min_alpha indexes
    for(auto n: leaves_to_root_iterator(tree, hg::leaves_it::exclude,hg::root_it::include)){
        if(min_alpha_marker[tree.child(0,n)] <= min_alpha_marker[tree.child(1,n)]){
            min_alpha_idx[n] = min_alpha_idx[tree.child(0, n)];
        }
        else{
            min_alpha_idx[n] = min_alpha_idx[tree.child(1, n)];
        }
    }
    closest_vertex[hg::root(tree)] = min_alpha_idx[hg::root(tree)];
    pass_edge_weight[hg::root(tree)] = altitudes[hg::root(tree)];
    //Checks from which node the min fminmax distance comes from.
    //Do that for each node on the tree
    for(auto n: root_to_leaves_iterator(tree, hg::leaves_it::include,hg::root_it::exclude)){
        auto dist_sib = min_alpha_marker[sibling[n]]*(1 - max_mu[sibling[n]] + altitudes[parents[n]]);
        if(fmm_distances[parents[n]] <= dist_sib){
            closest_vertex[n] = closest_vertex[parents[n]];
            fmm_distances[n] = fmm_distances[parents[n]];
            pass_edge_weight[n] = pass_edge_weight[parents[n]];
        }
        else{
            closest_vertex[n] = min_alpha_idx[sibling[n]];
            pass_edge_weight[n] = altitudes[parents[n]];
            fmm_distances[n] = dist_sib;
        }
    }
    //Checks if distance from leaf to all is not bigger than 
    //the distance from the leaves to themselves
    for(auto n: hg::leaves_iterator(tree)){
        float dist_to_itself = min_alpha_marker[n]*(1-marker[n]);
        if(dist_to_itself <= fmm_distances[n]){
            fmm_distances[n] = dist_to_itself;
            closest_vertex[n] = n;
            pass_edge_weight[n] = 0;
        }
    }
    return py::make_tuple(std::move(closest_vertex), std::move(pass_edge_weight));

}


// Example
auto example_function(const hg::tree & tree, const xt::pyarray<double> & altitudes){
    xt::pyarray<double> res = xt::zeros_like(altitudes);
    for(auto n: leaves_to_root_iterator(tree, hg::leaves_it::exclude)){
        double tmp = 1;
        for(auto c: hg::children_iterator(n, tree)){
            tmp *= res[c];
        }
        res[n] = altitudes[n] + tmp;
    }
    return res;
}

// Python Module and Docstrings
PYBIND11_MODULE(fminmax, m)
{
    xt::import_numpy();

    m.doc() = R"pbdoc(
        Calculation of the closest vertex for the Fminmax distance on a given tree.
    )pbdoc";
    m.def("example_function", example_function, "example_function");
    m.def("fminmax_idx", fminmax_idx, "Fminmax returining indexes of closest vertices and weight of the pass edges");
    
}
