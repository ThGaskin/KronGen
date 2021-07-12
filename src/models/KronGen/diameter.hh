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
  * \param diameter       The target diameter
  * \param log            The model logger
  */
template<typename Graph, typename Logger>
Graph create_first_Kronecker_factor(double& N_G,
                                    double& m_G,
                                    const std::string degree_distr,
                                    const double diameter,
                                    const Logger& log)
{
    log->info("Assembling diameter component, first factor ... ");

    // Create an empty graph
    Graph G{};

    // For scale-free base graphs, create a simple cycle of length 2*diameter
    // with degree variance 0
    if (degree_distr == "scale-free") {
        N_G = 2*diameter;
        m_G = 2;
        G = Utopia::Graph::create_regular_graph<Graph>(N_G, m_G, false);

        log->info("Result: cycle graph with {} vertices, mean degree {}.", N_G, m_G);
    }
    // For all others, create a chain of length diameter+1
    else {
        N_G = diameter + 1;
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
  * \param diameter       The target diameter
  * \param calculate_c    Whether or not to calculate the clustering coefficient
  * \param calculate_diam Whether or not to calculate the diameter
  * \param c_H            The clustering coefficient of the second factor
  * \param m_H            The mean_degree of the second factor
  * \param var_H          The degree distribution variance of the second factor
  * \param diam_H         The diameter of the second factor
  * \param discard_first_factor   Whether the first factor should be discarded
  * \param rng            The model rng
  * \param distr          Uniform distribution on (0, 1)
  * \param log            The model logger
  */
template<typename Graph, typename RNGType, typename Logger>
Graph create_second_Kronecker_factor(const double N,
                                     const double N_G,
                                     const double m,
                                     const double m_G,
                                     const std::string degree_distr,
                                     const double diameter,
                                     const bool calculate_c,
                                     const bool calculate_diam,
                                     double& c_H,
                                     double& m_H,
                                     double& var_H,
                                     double& diam_H,
                                     bool& discard_first_factor,
                                     RNGType& rng,
                                     std::uniform_real_distribution<double>& distr,
                                     const Logger& log)
{
    using namespace Utopia::Graph;

    log->info("Assembling diameter component, second factor ... ");

    // Create an empty graph
    Graph H{};
    double N_H = std::round(N/N_G);
    m_H = std::round((m+1)/(m_G+1)-1);

    // If no second Kronecker factor is possible: discard the first factor and
    // return
    if (m_H < 1) {
        discard_first_factor = true;
        return H;
    }

    // Get an estimate of the diameter of the second factor
    const auto estimated_diameter = std::round(Utils::diameter_estimation(N_H, m_H));

    // If the estimated diameter lower or equal to target diameter: create
    // second factor and combine.
    if (estimated_diameter <= diameter) {
        if (degree_distr == "scale-free") {
            if (m_H < 3) {
                discard_first_factor = true;
                log->info("Mean degree {} too small; discarding.", m_H);
                return H;
            }
            else {
                log->info("Result: scale-free graph with {} vertices, mean degree {}", N_H, m_H);
                H = create_BarabasiAlbert_graph<Graph>(N_H, m_H, false, rng);
            }
        }
        else {
            log->info("Result: random graph with {} vertices, mean degree {}", N_H, m_H);
            H = create_ErdosRenyi_graph<Graph>(N_H, m_H, false, false, rng);
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

    // Second Kronecker factor has diameter larger than estimate
    else {
        log->info("Diameter larger than estimate; attempting to further "
                  "break down second component ... ");
        // For small mean degree, no further breaking down of the mean degree is
        // possible
        if (m_H <= 2){
            discard_first_factor = true;
            log->info("Mean degree {} too small; discarding.", m_H);
            return H;
        }

        // Else: break down second factor into smaller units
        else {
            // Collect possible Kronecker factors of H
            std::vector<std::pair<std::size_t, double>> factors = {{N_H, m_H}};
            bool factors_too_large = true;

            while (factors_too_large) {
                factors_too_large = false;
                // Iterate over factors and break down any that are too large
                for (std::size_t i = 0; i < factors.size(); ++i) {
                    auto N_i = factors[i].first;
                    auto m_i = factors[i].second;
                    auto diam_i = std::round(Utils::diameter_estimation(N_i, m_i));

                    // Breaking cases:
                    if ( (N_i < 2)
                      or (m_i < 3 and degree_distr == "scale-free")
                      or ((N_i <= 2 or m_i <= 2) and (diam_i > diameter)))
                    {
                        log->info("Failed: N_i={}, m_i={}, diam_i={}; discarding.",
                                   N_i, m_i, diam_i);
                        discard_first_factor = true;
                        return H;
                    }

                    // Further break down factors
                    else if (diam_i > diameter) {
                        factors_too_large = true;
                        factors.erase(factors.begin()+i);

                        // Attempt to factor N_i
                        const auto new_factors = Utils::get_factors_N_m(N_i, m_i, rng);

                        if (new_factors.first.first == 0) {
                            discard_first_factor = true;
                            return H;
                        }

                        factors.push_back(new_factors.first);
                        factors.push_back(new_factors.second);

                        break;
                    }
                }
            }

            // If factors found, assemble H as Kronecker product of factors
            log->info("Factorisation successful; created {} factors",
                       factors.size());

            add_vertex(H);
            Utils::add_self_edges(H);

            for (const auto& f : factors) {
                Graph t{};
                if (degree_distr == "scale-free") {
                    t = create_BarabasiAlbert_graph<Graph>(f.first, f.second, false, rng);
                }
                else {
                      t = create_ErdosRenyi_graph<Graph>(f.first, f.second, false, false, rng);
                }

                if (calculate_c) {
                    const auto deg_stats = degree_statistics(t);
                    if (c_H == -1) {
                        c_H = global_clustering_coeff(t);
                        m_H = deg_stats.first;
                        var_H = deg_stats.second;
                    }
                    else {
                        const auto c_t = global_clustering_coeff(t);
                        const auto m_t = deg_stats.first;
                        const auto var_t = deg_stats.second;
                        c_H = Utils::Kronecker_clustering(c_H, c_t, m_H, m_t, var_H, var_t);
                        m_H = Utils::Kronecker_mean_degree(m_H, m_t);
                        var_H = Utils::Kronecker_degree_variance(m_H, m_t, var_H, var_t);
                    }
                }
                if (calculate_diam) {
                    diam_H = std::max(diam_H, Utopia::Models::NetworkAnalyser::diameter(t));
                }

                Utils::add_self_edges(t);

                H = Utils::Kronecker_product(H, t, rng, distr);

            }

            log->info("Done. Second diameter component properties: "
                      "N_H={}, m_H={}{}{}{}{}.",
                      N_H, m_H,
                      (calculate_c or calculate_diam ? ", " : ""),
                      (calculate_c ? "c_H="+to_string(c_H) : ""),
                      (calculate_c and calculate_diam ? ", " : ""),
                      (calculate_diam ? "diam_H="+to_string(diam_H) : ""));

            // Return the graph
            return H;
        }
    }
}

/// Create a graph with a given diameter
/**
  * \param K                The product graph
  * \param N                The target number of vertices
  * \param m                The target mean degree
  * \param diameter         The target diameter
  * \param degree_distr     The target degree distribution type
  * \param calculate_c      Whether or not to calculate the clustering coefficient
  * \param calculate_diam   Whether or not to calculate the diameter
  * \param rng              The model rng
  * \param distr            Uniform real distribution on (0, 1)
  * \param log              The model logger
  */
template<typename Graph, typename RNGType, typename Logger>
void create_diameter_graph (Graph& K,
                            const double N,
                            const double m,
                            const double c,
                            const double diameter,
                            const std::string degree_distr,
                            const bool calculate_c,
                            const bool calculate_diam,
                            RNGType& rng,
                            std::uniform_real_distribution<double>& distr,
                            const Logger& log)
{
    log->info("Assembling diameter component ... ");

    // Extreme case: mean degree less than 3
    if (m <= 3) {
        log->info("Mean degree < 3; returning star graph component with {} vertices, "
                  " mean degree {}", N, m);
        K = AuxGraphs::create_star_graph<Graph>(N, m, diameter, rng);
        K[0].state.diameter = diameter;
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
                                                   degree_distr, diameter,
                                                   log);
    const auto var_G = degree_statistics(G).second;

    // If clustering coefficient is specificed: return only first factor.
    if (c != -1) {
        K = G;
        K[0].state.mean_deg = m_G;
        K[0].state.var = var_G;
        K[0].state.clustering_global = c_G;
        K[0].state.diameter = diameter;

        log->info("Done: returning first diameter component.");

        return;
    }

    //... Second Kronecker factor ..............................................
    Graph H{};
    double c_H = -1, m_H = 0, var_H = 0, diam_H = 0;
    bool discard_first_factor = false;

    H = create_second_Kronecker_factor<Graph>(N, N_G, m, m_G,
                                                  degree_distr, diameter,
                                                  calculate_c, calculate_diam,
                                                  c_H, m_H, var_H, diam_H,
                                                  discard_first_factor,
                                                  rng, distr, log);

    // ... Create Kronecker graph and write properties .........................
    // Second factor successfully created
    if (not discard_first_factor) {

        K = Utils::Kronecker_product(G, H, rng, distr);

        K[0].state.mean_deg = Utils::Kronecker_mean_degree(m_G, m_H);
        K[0].state.var = Utils::Kronecker_degree_variance(m_G, m_H, var_G, var_H);
        if (calculate_c) {
            K[0].state.clustering_global =
                  Utils::Kronecker_clustering(c_G, c_H, m_G, m_H, var_G, var_H);
        }
        if (calculate_diam) {
            K[0].state.diameter = std::max(diameter, diam_H);
        }
    }

    // If H is empty, revert to star graph, dropping first factor
    else {
        log->info("Discarding first factor and reverting to star graph with "
                  "N={}, m={}.", N, m);

        K = AuxGraphs::create_star_graph<Graph>(N, m, diameter, rng);
        K[0].state.mean_deg = m;
        K[0].state.var = degree_statistics(K).second;
        if (calculate_c) {
            K[0].state.clustering_global = global_clustering_coeff(K);
        }
        if (calculate_diam) {
            K[0].state.diameter = diameter;
        }
    }
}

}
#endif
