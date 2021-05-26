#ifndef UTOPIA_MODELS_KRONGEN_GRAPHCREATION
#define UTOPIA_MODELS_KRONGEN_GRAPHCREATION

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/random.hpp>
#include <boost/property_map/dynamic_property_map.hpp>

#include "utopia/data_io/cfg_utils.hh"
#include "utopia/data_io/filesystem.hh"
#include "utopia/data_io/graph_load.hh"
#include "utopia/core/types.hh"
#include "utopia/core/graph/iterator.hh"
#include "utopia/core/graph/creation.hh"

namespace Utopia::Models::KronGen::GraphCreation {

using namespace Utopia;
using namespace boost;

// Kronecker product of graphs
// Graphs must have a self-loop on every node
template<typename Graph>
Graph Kronecker_product(Graph& K, Graph& G) {
    Graph P{boost::num_vertices(K)*boost::num_vertices(G)};
    const std::size_t M = boost::num_vertices(G);
    for(const auto k : range<IterateOver::edges>(K)) {
        for(const auto g : range<IterateOver::edges>(G)) {
            const auto s_1 = source(k, K);
            const auto s_2 = source(g, G);
            const auto t_1 = target(k, K);
            const auto t_2 = target(g, G);

            auto i = s_1*M+s_2;
            auto j = t_1*M+t_2;
            if ((not edge(i, j, P).second) && (not edge(j, i, P).second)){
                add_edge(i, j, P);
            }

            i = s_1*M+t_2;
            j = t_1*M+s_2;
            if ((not edge(i, j, P).second) && (not edge(j, i, P).second)){
                add_edge(i, j, P);
            }
        }
    }
    return P;
}

template<typename Graph, typename RNGType>
Graph create_Kronecker_graph(const Config& cfg,
                             RNGType& rng)
{
    // Get the graph generating model
    const std::string model = get_as<std::string>("model", cfg);

    if (model != "Kronecker") {
        return Utopia::Graph::create_graph<Graph>(cfg, rng);
    }

    // Create empty graph and add one vertex with a self-loop
    Graph g{};
    auto v = add_vertex(g);
    add_edge(v, v, g);

    // Generate Graph
    for (const auto& model_map : cfg["Kronecker"]){
        const auto& model_cfg = model_map.second;
        Graph h = Utopia::Graph::create_graph<Graph>(model_cfg, rng);
        for (const auto v : range<IterateOver::vertices>(h)) {
            add_edge(v, v, h);
        }
        g = Kronecker_product(g, h);
    }

    // Remove self-loops
    for (const auto v : range<IterateOver::vertices>(g)) {
        remove_edge(v, v, g);
    }

    // Return the graph
    return g;
}
}
#endif
