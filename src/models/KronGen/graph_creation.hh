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

#include "../NetworkAnalyser/graph_metrics.hh"
#include "utils.hh"

namespace Utopia::Models::KronGen::GraphCreation {

using namespace boost;


/// Create a graph with zero clustering, a given number of vertices, a
/// given mean degree, and zero degree variance.
/**
  * \param N    Number of vertices. Must be at least 2*k, and even.
  * \param k    Mean degree of graph. Must be at least 2.
  * \return g   The graph
*/
template<typename Graph>
Graph create_zero_c_graph(const std::size_t N, const std::size_t k)
{

    if (N < 2*k) {
        throw std::invalid_argument("N must be at least equal to 2k!");
    }

    else if (N%2) {
        throw std::invalid_argument("N must be even!");
    }

    if (k < 2) {
        throw std::invalid_argument("k must be at least 2!");
    }

    Graph g{N};

    for (std::size_t v = 0; v < N/2; ++v) {
        const auto v_2 = v + N/2 - k;
        for (std::size_t i = 0; i < k; ++i) {
            auto nb = (v_2 + i)%(N/2) + N/2;
            add_edge(v, nb, g);
        }
    }

    return g;
}

// Kronecker product of graphs
// Graphs must have a self-loop on every node
// Todo: Test me
template<typename Graph>
Graph Kronecker_product(Graph& K, Graph& G) {
    using namespace Utopia;

    const std::size_t M = num_vertices(G);
    Graph P{num_vertices(K)*M};

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
Graph create_Kronecker_graph(const Config& cfg, RNGType& rng)
{

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

template<typename Graph, typename RNGType>
Graph create_KronGen_graph(const Config& cfg, RNGType& rng)
{

    std::uniform_real_distribution<double> distr(0, 1);

    // Create an empty graph
    Graph g{};
    auto v = add_vertex(g);
    add_edge(v, v, g);

    const double c = get_as<double>("clustering_coeff", cfg["KronGen"]);
    std::size_t N;
    std::size_t m;

    // Fixme: do this properly
    // For low clustering: low k/N, but with sufficiently large m
    if (c < 0.2) {
        N = std::round(get_as<std::size_t>("num_vertices", cfg)*0.8);
        m = 8;
    }
    // For high clustering: high k/N.
    // for very high clustering: reduce k slightly
    else {
        N = std::round(get_as<std::size_t>("num_vertices", cfg)*0.2);
        if (c > 0.8) {
            m = std::round(N-3);
        }
        else {
            m = std::round(N-1);
        }
    }

    Graph r = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N, m, false, false, rng);
    double c_r = Utopia::Models::NetworkAnalyser::global_clustering_coeff(r);
    auto deg = Utopia::Models::NetworkAnalyser::degree_statistics(r);
    double mean_deg = deg.first;
    double var = deg.second;

    for (const auto v : range<IterateOver::vertices>(r)) {
        add_edge(v, v, r);
    }
    if (c_r <= c) {
        const std::size_t n_2 = Utils::get_mean_deg_c(false, c_r, c, mean_deg, var);
        Graph k = Utopia::Graph::create_complete_graph<Graph>(n_2+1);
        for (const auto v : range<IterateOver::vertices>(k)) {
            add_edge(v, v, k);
        }
        r = Kronecker_product(k, r);
        c_r = Utils::Kronecker_clustering(c_r, 1, mean_deg, n_2, var, 0);
    }
    else {
        const std::size_t n_2 = Utils::get_mean_deg_c(true, c_r, c, mean_deg, var);
        Graph k = create_zero_c_graph<Graph>(2*n_2, n_2);
        for (const auto v : range<IterateOver::vertices>(k)) {
            add_edge(v, v, k);
        }
        r = Kronecker_product(k, r);
        c_r = Utils::Kronecker_clustering(c_r, 0, mean_deg, n_2, var, 0);
    }
    for (const auto v : range<IterateOver::vertices>(r)) {
        remove_edge(v, v, r);
    }

    return r;

}

// Custom create_graph function
template<typename Graph, typename RNGType>
Graph create_graph(const Config& cfg, RNGType& rng)
{
  // Get the graph generating model
  const std::string model = get_as<std::string>("model", cfg);

  std::uniform_real_distribution<double> distr(0, 1);
  if (model == "Kronecker") {
      return create_Kronecker_graph<Graph>(cfg, rng);
  }
  else if (model == "KronGen") {
      return create_KronGen_graph<Graph>(cfg, rng);
  }

  else {
      return Utopia::Graph::create_graph<Graph>(cfg, rng);
  }
}
}
#endif
