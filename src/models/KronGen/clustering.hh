#ifndef UTOPIA_MODELS_KRONGEN_CLUSTERING
#define UTOPIA_MODELS_KRONGEN_CLUSTERING

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/random.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <spdlog/spdlog.h>

#include "utopia/core/types.hh"
#include "utopia/core/graph/iterator.hh"
#include "utopia/core/graph/creation.hh"
#include "utopia/data_io/cfg_utils.hh"
#include "utopia/data_io/filesystem.hh"
#include "utopia/data_io/graph_load.hh"

#include "aux_graphs.hh"
#include "utils.hh"

#include "../NetworkAnalyser/graph_metrics.hh"

namespace Utopia::Models::KronGen::Clustering {

using namespace boost;
using namespace Utopia::Models::KronGen;
using namespace Utopia::Models::NetworkAnalyser;

/// Adjust the size of the initial network to the clustering coefficient
// To do: Do this properly
template<typename Logger>
void adjust_N_m_to_c(double& N, double& m, const double c, const Logger& log)
{
    // If clustering coefficient is not specified, do nothing
    if (c == -1) {
        return;
    }

    // For low clustering: low k/N, but with sufficiently large m
    else if (c < 0.2) {
        N = std::round(N*0.8);
        m = 8;
    }

    // For high clustering: high k/N.
    // for very high clustering: reduce k slightly
    else {
        N = std::round(N*0.2);
        if (c > 0.8) {
            m = std::round(N-3);
        }
        else {
            m = std::round(N-1);
        }
    }

    log->info("Adjusted N to {}, m to {}; (c={}).", N, m, c);

}

/// Create a graph with a given clustering
template<typename Graph, typename RNGType, typename Logger>
void create_clustering_graph(Graph& K,
                             double N,
                             double m,
                             const double c,
                             const double diameter,
                             const std::string degree_distr,
                             const bool calculate_c,
                             const bool calculate_diam,
                             RNGType& rng,
                             std::uniform_real_distribution<double>& distr,
                             const Logger& log)
{

    using vertices_size_type = typename graph_traits<Graph>::vertices_size_type;

    log->info("Assembling clustering component ... ");

    // ... Get first Kronecker factor ..........................................
    Graph G{};
    double c_G, diam_G, m_G, var_G;

    if (diameter > 0) {
        log->info("First factor from previous assembly.");
        G = K;  // already has self-edges on every vertex
        c_G = K[0].state.clustering_global;
        diam_G = K[0].state.diameter;
        m_G = K[0].state.mean_deg;
        var_G = K[0].state.var;
    }
    else {
        log->info("Assembling clustering component, first factor ...");
        adjust_N_m_to_c(N, m, c, log);
        auto N_G = N;
        m_G = m;
        if (degree_distr == "scale-free") {
            G = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N_G, m_G, false, rng);
            log->info("Result: scale-free graph with {} vertices, mean degree {}.", N_G, m_G);
        }
        else {
            G = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N_G, m_G, false, false, rng);
            log->info("Result: random graph with {} vertices, mean degree {}.", N_G, m_G);
        }

        c_G = global_clustering_coeff(G);
        const auto deg_stats = degree_statistics(G);
        m_G = deg_stats.first;
        var_G = deg_stats.second;
        if (calculate_diam) {
            const auto s = fourSweep<vertices_size_type>(G);
            diam_G = iFUB(s.first, s.second, 0, G);
        }

        log->info("First factor properties: c_G={}{}.", c_G,
                  (calculate_diam ? ", diam_G="+to_string(diam_G) : ""));

        Utils::add_self_edges(G);
    }

    // ... Get second Kronecker factor .........................................
    log->info("Assembling clustering component, second factor ...");
    double N_H, c_H, diam_H, m_H, var_H;
    m_H = std::round(Utils::get_mean_deg_c(c_G, c, m_G, var_G));

    // FIXME adjust this .......................................................
    if (c_G <= c) {
        N_H = m_H +1;
    }
    else {
        N_H = std::max(6., 2*m_H);
    }
    Graph H = (c_G <= c)
      ? Utopia::Graph::create_complete_graph<Graph>(N_H)
      : AuxGraphs::create_zero_c_graph<Graph>(N_H, m_H);

    log->info("Result: {} graph with N_H={}, m_H={}", (c_G <= c ? "complete" : "zero_c"), N_H, m_H);

    // .........................................................................

    // Calculate second factor properties
    c_H = (c_G <= c) ? 1 : 0;
    var_H = 0;
    if (calculate_diam) {
        const auto s = fourSweep<vertices_size_type>(H);
        diam_H = iFUB(s.first, s.second, 0, H);
        log->info("Second factor diameter: {}", diam_H);
    }

    // ... Create Kronecker graph and write properties .........................
    Utils::add_self_edges(H);

    K = Utils::Kronecker_product(G, H, rng, distr);

    if (calculate_c) {
        K[0].state.clustering_global = Utils::Kronecker_clustering(c_G, c_H,
                                                                   m_G, m_H,
                                                                   var_G, var_H);
        log->info("Kronecker product clustering_coeff: {}", K[0].state.clustering_global);
    }

    if (calculate_diam) {
        K[0].state.diameter = std::max(diam_G, diam_H);
        log->info("Kronecker product diameter: {}", K[0].state.diameter);
    }

}
}
#endif
