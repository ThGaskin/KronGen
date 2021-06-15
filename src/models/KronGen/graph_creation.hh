#ifndef UTOPIA_MODELS_KRONGEN_GRAPHCREATION
#define UTOPIA_MODELS_KRONGEN_GRAPHCREATION

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/random.hpp>
#include <boost/property_map/dynamic_property_map.hpp>

#include "utopia/core/types.hh"
#include "utopia/core/graph/iterator.hh"
#include "utopia/core/graph/creation.hh"
#include "utopia/data_io/cfg_utils.hh"
#include "utopia/data_io/filesystem.hh"
#include "utopia/data_io/graph_load.hh"

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

/// Create a 1-d lattice chain of length N. The diameter of the chain is N-1
/**
  * \param N   Length of chain
*/
template<typename Graph>
Graph create_chain_graph(const std::size_t N)
{
    Graph g = Utopia::Graph::create_regular_graph<Graph>(N, 2, false);
    remove_edge(0, 1, g);

    return g;
}

/// Returns the Kronecker product of a list of graphs
template<typename Graph, typename RNGType>
Graph create_Kronecker_graph(const Config& cfg,
                             RNGType& rng,
                             const Config& analysis_cfg = YAML::Node(YAML::NodeType::Map))
{
    // Create empty graph and add one vertex with a self-loop
    Graph g{};
    auto v = add_vertex(g);
    add_edge(v, v, g);

    // Data containers for graph analysis properties
    double c_global = -1;
    double diam = -1;

    double mean_deg;
    double variance;
    bool calculate_c = false;
    bool calculate_diam = false;
    try {
       calculate_c = (get_as<bool>("clustering_global", analysis_cfg)
                    && get_as<bool>("enabled", analysis_cfg));
       calculate_diam = (get_as<bool>("diameter", analysis_cfg)
                    && get_as<bool>("enabled", analysis_cfg));
    }
    catch (YAML::InvalidNode) {}
    catch (Utopia::KeyError) {}
    bool first_run = true;

    // Generate Graph
    for (const auto& model_map : cfg["Kronecker"]){
        const auto& model_cfg = model_map.second;
        const auto& model = get_as<std::string>("model", model_cfg);
        Graph h;

        if (model == "zero_c") {
            h = create_zero_c_graph<Graph>(
                get_as<std::size_t>("num_vertices", model_cfg),
                get_as<std::size_t>("mean_degree", model_cfg)
            );
        }
        else if (model == "chain") {
            h = create_chain_graph<Graph>(
                get_as<std::size_t>("num_vertices", model_cfg)
            );
        }
        else {
            h = Utopia::Graph::create_graph<Graph>(model_cfg, rng);
        }

        // ... Calulate properties: stored in first vertex .....................
        // Global clustering coefficient
        Utils::calculate_properties (g,
                                     h,
                                     first_run,
                                     calculate_c,
                                     calculate_diam,
                                     c_global,
                                     diam,
                                     mean_deg,
                                     variance);
        // .....................................................................

        for (const auto v : range<IterateOver::vertices>(h)) {
            add_edge(v, v, h);
        }

        g = Utils::Kronecker_product(g, h);
    }

    // Remove self-loops
    for (const auto v : range<IterateOver::vertices>(g)) {
        remove_edge(v, v, g);
    }

    // Write data
    if (calculate_c) {
        g[0].state.clustering_global = c_global;
    }
    if(calculate_diam) {
        g[0].state.diameter = diam;
    }

    // Return the graph
    return g;
}

// Creates a graph from a list of topological properties
template<typename Graph, typename RNGType>
Graph create_KronGen_graph(const Config& cfg,
                           RNGType& rng,
                           const Config& analysis_cfg = YAML::Node(YAML::NodeType::Map))
{

    using vertices_size_type = typename boost::graph_traits<Graph>::vertices_size_type;

    std::uniform_real_distribution<double> distr(0, 1);

    // Create an empty graph
    Graph g{};
    auto v = add_vertex(g);
    add_edge(v, v, g);

    const double c = get_as<double>("clustering_coeff", cfg["KronGen"]);
    std::size_t N;
    std::size_t m;

    // Fixme: do this properly..................................................
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
    // .........................................................................

    Graph r = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N, m, false, false, rng);
    double c_r = Utopia::Models::NetworkAnalyser::global_clustering_coeff(r);
    auto deg = Utopia::Models::NetworkAnalyser::degree_statistics(r);
    double mean_deg = deg.first;
    double var = deg.second;
    double diam = -1;
    const bool calculate_diam = (analysis_cfg != YAML::Node(YAML::NodeType::Map))
      ? (get_as<bool>("diameter", analysis_cfg)
         && get_as<bool>("enabled", analysis_cfg))
      : false;
    if (calculate_diam) {
        const auto starting_point = Utopia::Models::NetworkAnalyser::fourSweep<vertices_size_type>(r);
        diam = Utopia::Models::NetworkAnalyser::iFUB(starting_point.first, starting_point.second, 0, r);
    }

    for (const auto v : range<IterateOver::vertices>(r)) {
        add_edge(v, v, r);
    }

    const std::size_t n_2 = std::round(Utils::get_mean_deg_c(c_r, c, mean_deg, var));
    Graph k = (c_r <= c) ? Utopia::Graph::create_complete_graph<Graph>(n_2+1) : create_zero_c_graph<Graph>(2*n_2, n_2);
    if (calculate_diam) {
        const auto starting_point = Utopia::Models::NetworkAnalyser::fourSweep<vertices_size_type>(k);
        diam = std::max(static_cast<double>(Utopia::Models::NetworkAnalyser::iFUB(starting_point.first, starting_point.second, 0, k)), diam);
    }
    for (const auto v : range<IterateOver::vertices>(k)) {
        add_edge(v, v, k);
    }
    r = Utils::Kronecker_product(k, r);
    const double c_K = (c_r <= c) ? 1 : 0;
    c_r = Utils::Kronecker_clustering(c_r, c_K, mean_deg, n_2, var, 0);
    r[0].state.clustering_global = c_r;
    r[0].state.diameter = diam;

    for (const auto v : range<IterateOver::vertices>(r)) {
        remove_edge(v, v, r);
    }

    return r;
}

// Custom create_graph function
template<typename Graph, typename RNGType>
Graph create_graph(const Config& cfg, RNGType& rng, const bool includes_analysis_cfg = false)
{
  const Config graph_cfg = includes_analysis_cfg
      ? get_as<Config>("create_graph", cfg)
      : cfg;
  Config nw_cfg;

  try {
      nw_cfg = get_as<Config>("graph_analysis", cfg["NetworkAnalyser"]);
  }
  catch (YAML::InvalidNode){
      nw_cfg = YAML::Node(YAML::NodeType::Map);
  }
  catch (Utopia::KeyError){
      nw_cfg = YAML::Node(YAML::NodeType::Map);
  }
  // Get the graph generating model
  const std::string model = get_as<std::string>("model", graph_cfg);

  std::uniform_real_distribution<double> distr(0, 1);
  if (model == "Kronecker") {
      return create_Kronecker_graph<Graph>(graph_cfg, rng, nw_cfg);
  }
  else if (model == "KronGen") {
      return create_KronGen_graph<Graph>(graph_cfg, rng, nw_cfg);
  }

  else if (model == "zero_c") {
      return create_zero_c_graph<Graph>(
          get_as<std::size_t>("num_vertices", graph_cfg),
          get_as<std::size_t>("mean_degree", graph_cfg)
      );
  }
  else if (model == "chain") {
      return create_chain_graph<Graph>(get_as<std::size_t>("num_vertices", graph_cfg));
  }

  else {
      return Utopia::Graph::create_graph<Graph>(graph_cfg, rng);
  }
}
}
#endif
