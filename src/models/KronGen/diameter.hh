#ifndef UTOPIA_MODELS_KRONGEN_DIAMETER
#define UTOPIA_MODELS_KRONGEN_DIAMETER

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/random.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <spdlog/spdlog.h>

#include "aux_graphs.hh"
#include "clustering.hh"
#include "utils.hh"

#include "../NetworkAnalyser/graph_metrics.hh"

namespace Utopia::Models::KronGen::Diameter {

using namespace boost;
using namespace Utopia::Models::KronGen;
using namespace Utopia::Models::NetworkAnalyser;

/// ... Methods used in creating a graph with a given diameter .................

/// Create the first Kronecker factor. Returns a graph with self-edges on every
/// node.
/**
  * \param N_G            The number of vertices of the first factor
  * \param m_G            The mean degree of the first factor
  * \param degree_distr   The target degree distribution
  * \param d              The target diameter
  * \param log            The model logger
  */
template<typename Graph, typename Logger>
Graph create_first_Kronecker_factor(double& N_G,
                                    double& m_G,
                                    const std::string degree_distr,
                                    const double d,
                                    const Logger& log)
{
    log->info("Assembling diameter component, first factor ... ");

    // Create an empty graph
    Graph G{};

    // For scale-free base graphs, create a simple cycle of length 2*diameter
    // with degree variance 0
    if (degree_distr == "scale-free") {
        N_G = 2*d;
        m_G = 2;
        G = Utopia::Graph::create_regular_graph<Graph>(N_G, m_G, false);

        log->info("Result: cycle graph with {} vertices, mean degree {}.", N_G, m_G);
    }
    // For all others, create a chain of length diameter+1
    else {
        N_G = d + 1;
        m_G = Utils::mean_degree_chain_graph(N_G);
        G = AuxGraphs::create_chain_graph<Graph>(N_G);

        log->info("Result: chain graph with {} vertices, mean degree {}.", N_G, m_G);
    }

    Utils::add_self_edges(G);

    // Return the graph
    return G;
}

/// Create the second Kronecker factor. Returns a graph with a self-edge on
/// every node.
/**
  * \param N              The number of vertices of the final graph
  * \param N_G            The number of vertices of the first factor
  * \param m              The mean degree of the final graph
  * \param m_G            The mean degree of the first factor
  * \param degree_distr   The target degree distribution
  * \param d              The target diameter
  * \param calculate_c    Whether or not to calculate the clustering coefficient
  * \param calculate_diam Whether or not to calculate the diameter
  * \param tolerance      The tolerance value for the grid search
  * \param c_H            The clustering coefficient of the second factor
  * \param m_H            The mean_degree of the second factor
  * \param var_H          The degree distribution variance of the second factor
  * \param diam_H         The diameter of the second factor
  * \param discard_first_factor   Whether the first factor should be discarded
  * \param rng            The model rng
  * \param log            The model logger
  */
template<typename Graph, typename RNGType, typename Logger>
Graph create_second_Kronecker_factor(const double N,
                                     const double N_G,
                                     const double m,
                                     const double m_G,
                                     const std::string degree_distr,
                                     const double d,
                                     const bool calculate_c,
                                     const bool calculate_diam,
                                     const double tolerance,
                                     double& c_H,
                                     double& m_H,
                                     double& var_H,
                                     double& diam_H,
                                     bool& discard_first_factor,
                                     RNGType& rng,
                                     const Logger& log)
{
    using namespace Utopia::Graph;

    log->info("Assembling diameter component, second factor ... ");

    // Create an empty graph
    Graph H{};
    double N_H = std::round(N/N_G);
    m_H = std::round(Utils::Kronecker_mean_degree_inv(m, m_G));

    log->debug("Required: N_H = {}, m_H = {}.", N_H, m_H);

    // If no second Kronecker factor is possible: discard the first factor and
    // return
    if (m_H < 1 or (degree_distr == "scale-free" and m_H < 3)) {
        discard_first_factor = true;
        return H;
    }

    // Get an estimate of the diameter of the second factor
    const auto estimated_diameter = std::round(Utils::diameter_estimation(N_H, m_H));
    log->debug("Estimated diameter = {}", estimated_diameter);

    // If the estimated diameter lower or equal to target diameter: create
    // second factor and combine.
    // To do: this can all be integrated into one grid search
    if (estimated_diameter <= d) {
        if (degree_distr == "scale-free") {
            log->info("Result: scale-free graph with {} vertices, mean degree {}", N_H, m_H);
            H = create_BarabasiAlbert_graph<Graph>(N_H, m_H, false, rng);
        }
        else {
            log->info("Result: random graph with {} vertices, mean degree {}", N_H, m_H);
            H = create_ErdosRenyi_graph<Graph>(N_H, m_H, false, false, rng);
        }
    }

    // Second Kronecker factor has diameter larger than estimate
    else {
        log->info("Assembling second component via grid search; component "
                  "target size: N_H={}, m_H={}", N_H, m_H);

        // Get factors of N_H, m_H
        auto N_factors = Utils::closest_N_factors(N_H);
        for (const auto& nn : N_factors) {
            N_factors.push_back({nn.second, nn.first});
        }
        auto m_factors = Utils::closest_mean_deg_factors(m_H);
        for (const auto& mm : m_factors) {
            m_factors.push_back({mm.second, mm.first});
        }

        double diam_H = estimated_diameter;

        // Base error
        double error = Clustering::err_func({
            Clustering::rel_err(N_G*N_H, N),
            Clustering::rel_err((m_G+1)*(m_H+1), m+1),
            Clustering::rel_err(std::max(diam_H, d), d)});

        double N_1=N_H, m_1=m_H, diam_1=diam_H;
        double N_2=1, m_2=0, diam_2=0;

        // ... Grid search .........................................................
        if (error > tolerance) {
            log->debug("Commencing grid search ...");

            // Current error value
            double current_err;

            for (const auto& m_fac: m_factors) {

                for (const auto& n_fac : N_factors) {

                    if ((m_fac.first >= n_fac.first)
                       or (m_fac.second >= n_fac.second)
                       or (degree_distr == "scale-free" and m_fac.first < 3)
                       or (degree_distr == "scale-free" and m_fac.second < 3)) {
                          continue;
                    }

                    log->debug("Checking case N_H={}x{}, m_H={}x{} ... ",
                          n_fac.first, n_fac.second, m_fac.first, m_fac.second);

                    diam_1 = Utils::diameter_estimation(n_fac.first, m_fac.first);
                    diam_2 = Utils::diameter_estimation(n_fac.second, m_fac.second);

                    log->debug("Estimated diameters: {}, {}", diam_1, diam_2);

                    // Calculate error
                    current_err = Clustering::err_func(
                        {Clustering::rel_err(n_fac.first*n_fac.second*N_G, N),
                         Clustering::rel_err((m_fac.first+1)*(m_fac.second+1)*(m_G+1), m+1),
                         Clustering::rel_err(std::max({diam_1, diam_2, d}), d)
                        });

                    // If the current graph reduces the error: set graph factor
                    // values to current state
                    if (current_err < error) {
                        N_1 = n_fac.first;
                        N_2 = n_fac.second;
                        m_1 = m_fac.first;
                        m_2 = m_fac.second;
                        diam_H = std::max(diam_1, diam_2);
                        error = current_err;
                        log->debug("Improvement: N_H={}x{}, m_H={}x{}, diam_H={}",
                          N_1, N_2, m_1, m_2, diam_H);
                    }
                    if (error < tolerance) {break;}
                }
                if (error < tolerance) {break;}
            }

            log->debug("Grid search complete. Preliminary results: N_H={}x{}, "
                       "m_H={}x{}, diam_H = {}", N_1, N_2, m_1, m_2, diam_H);
        }

        log->info("Generating graph H ...");

        // Generate the components G and H based on the values obtained from the
        // grid search
        if (N_1 == N_H) {
            if (degree_distr == "scale-free") {
                H = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N_H, m_H,
                                                                      false, rng);
            }
            else {
                H = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N_H, m_H,
                                                             false, false, rng);
            }
        }
        else {
            Graph T1, T2;
            if (degree_distr == "scale-free") {
                T1 = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N_1, m_1,
                                                                      false, rng);
                T2 = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N_2, m_2,
                                                                      false, rng);
            }
            else {
                T1 = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N_1, m_1,
                                                             false, false, rng);
                T2 = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N_2, m_2,
                                                             false, false, rng);
            }
            Utils::add_self_edges(T1);
            Utils::add_self_edges(T2);
            H = Utils::Kronecker_product(T1, T2);
            Utils::remove_self_edges(H);
        }
    }
    const auto deg_stats = degree_statistics(H);
    m_H = deg_stats.first;
    var_H = deg_stats.second;

    if (calculate_c) {
        c_H = global_clustering_coeff(H);
    }

    if (calculate_diam) {
        diam_H = Utopia::Models::NetworkAnalyser::diameter(H);
    }
    log->info("Second diameter component properties: {}{}{}",
              (calculate_c ? "c_H="+to_string(c_H) : ""),
              (calculate_c and calculate_diam ? ", " : ""),
              (calculate_diam ? "diam_H="+to_string(diam_H) : ""));

    Utils::add_self_edges(H);

    return H;
}

/// Create a graph with a given diameter
/**
  * \param K                The product graph
  * \param N                The target number of vertices
  * \param m                The target mean degree
  * \param d                The target diameter
  * \param degree_distr     The target degree distribution type
  * \param calculate_c      Whether or not to calculate the clustering coefficient
  * \param calculate_diam   Whether or not to calculate the diameter
  * \param rng              The model rng
  * \param log              The model logger
  */
template<typename Graph, typename RNGType, typename Logger>
void create_diameter_graph (Graph& K,
                            const double N,
                            const double m,
                            const double c,
                            const double d,
                            const std::string degree_distr,
                            const bool calculate_c,
                            const bool calculate_diam,
                            const double tolerance,
                            RNGType& rng,
                            const Logger& log)
{
    log->info("Assembling diameter component ... ");

    // Extreme case: mean degree less than 3
    if (m <= 3.0) {
        log->info("Mean degree < 3; returning star graph component with {} vertices, "
                  " mean degree {}", N, m);
        K = AuxGraphs::create_star_graph<Graph>(N, m, d, rng);
        K[0].state.diameter = d;
        if (calculate_c) {
            K[0].state.clustering_global = global_clustering_coeff(K);
        }

        return;
    }

    //... First Kronecker factor ...............................................
    const auto c_G = 0;
    double N_G;
    double m_G;
    Graph G = create_first_Kronecker_factor<Graph>(N_G, m_G,
                                                   degree_distr, d,
                                                   log);
    const auto var_G = degree_statistics(G).second;

    // If clustering coefficient is specificed: return only first factor.
    if (c != -1) {
        K = G;
        K[0].state.mean_deg = m_G;
        K[0].state.var = var_G;
        K[0].state.clustering_global = c_G;
        K[0].state.diameter = d;

        log->info("Done: returning first diameter component.");

        return;
    }

    //... Second Kronecker factor ..............................................
    Graph H{};
    double c_H = -1, m_H = 0, var_H = 0, diam_H = 0;
    bool discard_first_factor = false;

    H = create_second_Kronecker_factor<Graph>(N, N_G, m, m_G,
                                                  degree_distr, d,
                                                  calculate_c, calculate_diam, tolerance,
                                                  c_H, m_H, var_H, diam_H,
                                                  discard_first_factor,
                                                  rng, log);

    // ... Create Kronecker graph and write properties .........................
    // Second factor successfully created
    if (not discard_first_factor) {

        K = Utils::Kronecker_product(G, H);

        K[0].state.mean_deg = Utils::Kronecker_mean_degree(m_G, m_H);
        K[0].state.var = Utils::Kronecker_degree_variance(m_G, m_H, var_G, var_H);
        if (calculate_c) {
            K[0].state.clustering_global =
                  Utils::Kronecker_clustering(c_G, c_H, m_G, m_H, var_G, var_H);
        }
        if (calculate_diam) {
            K[0].state.diameter = std::max(d, diam_H);
        }
    }

    // If H is empty, revert to star graph, dropping first factor
    else {
        log->info("Discarding first factor and reverting to star graph with "
                  "N={}, m={}.", N, m);

        K = AuxGraphs::create_star_graph<Graph>(N, m, d, rng);
        K[0].state.mean_deg = m;
        K[0].state.var = degree_statistics(K).second;
        if (calculate_c) {
            K[0].state.clustering_global = global_clustering_coeff(K);
        }
        if (calculate_diam) {
            K[0].state.diameter = d;
        }
    }
}

}
#endif
