#ifndef UTOPIA_MODELS_KRONGEN_CLUSTERING
#define UTOPIA_MODELS_KRONGEN_CLUSTERING

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

#include "aux_graphs.hh"
#include "utils.hh"

#include "../NetworkAnalyser/graph_metrics.hh"

namespace Utopia::Models::KronGen::Clustering {

using namespace boost;
using namespace Utopia::Models::KronGen;
using namespace Utopia::Models::NetworkAnalyser;

template<typename Graph, typename RNGType>
void create_clustering_graph(Graph& K,
                             double N,
                             double m,
                             const double c,
                             const double diameter,
                             const std::string degree_distr,
                             const bool calculate_c,
                             const bool calculate_diam,
                             RNGType& rng)
{

    using vertices_size_type = typename boost::graph_traits<Graph>::vertices_size_type;

    // Trash: unused parameters
    if (degree_distr == "fixme") {
        auto trash = diameter;
        trash += 1;
    }

    // Fixme: do this properly..............................................
    // For low clustering: low k/N, but with sufficiently large m
    if (c < 0.2) {
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
    // .........................................................................

    Graph r{};
    r = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N, m, false, false, rng);
    double c_r = global_clustering_coeff(r);
    auto deg = degree_statistics(r);
    double mean_deg = deg.first;
    double var = deg.second;
    double diam = 0;

    if (calculate_diam) {
        const auto starting_point = fourSweep<vertices_size_type>(r);
        diam = iFUB(starting_point.first, starting_point.second, 0, r);
    }

    Utils::add_self_edges(r);

    const std::size_t n_2 = std::round(Utils::get_mean_deg_c(c_r, c, mean_deg, var));

    // FIXME adjust this .......................................................
    Graph k = (c_r <= c)
      ? Utopia::Graph::create_complete_graph<Graph>(n_2+1)
      : AuxGraphs::create_zero_c_graph<Graph>(std::max(static_cast<std::size_t>(6), 2*n_2), n_2);
    // .........................................................................

    if (calculate_diam) {
        const auto starting_point =
            Utopia::Models::NetworkAnalyser::fourSweep<vertices_size_type>(k);
        diam = std::max(
                static_cast<double>(
                  Utopia::Models::NetworkAnalyser::iFUB(starting_point.first,
                                                        starting_point.second,
                                                        0,
                                                        k)),
                diam
        );
    }
    Utils::add_self_edges(k);
    K = Utils::Kronecker_product(k, r);
    const double c_K = (c_r <= c) ? 1 : 0;
    if (calculate_c) {
        c_r = Utils::Kronecker_clustering(c_r, c_K, mean_deg, n_2, var, 0);
        K[0].state.clustering_global = c_r;
    }

    K[0].state.diameter = std::max(K[0].state.diameter, diam);

    Utils::remove_self_edges(K);

}
}
#endif
