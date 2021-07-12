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

// ... Grid search utility functions ...........................................

/// Relative error of variable x to y
double rel_err(const double x, const double y) {

    return abs(1.-x/y);
}

/// Square mean of a list of errors
double err_func(const std::vector<double> err_terms) {
    double err = 0;
    for (const auto& x : err_terms) {
        err += pow(x, 2);
    }

    return sqrt(err);
}

/// Diameter error
double diameter_error(const double diam_1, const double diam_target) {
    const double err = (diam_1 > diam_target and diam_target > 0)
      ? rel_err(diam_1, diam_target)
      : 0;

    return err;
}

/// Create a graph with a given clustering
/**
  *\ param K                The graph to be returned, passed down from previous
                            proccesses
  *\ param N                The target number of vertices
  *\ param m                The target mean degree
  *\ param c                The target clustering coefficient
  * \param diameter         The target diameter
  * \param degree_distr     The target degree distribution
  * \param calculate_c      Whether or not to calculate the clustering coefficient
  * \param calculate_diam   Whether or not to calculate the diameter
  * \param tolerance        Tolerance in adjusting the graph properties
  * \param rng              The model rng
  * \param distr            Uniform real distribution
  * \param log              The model loggger
 */
template<typename Graph, typename RNGType, typename Logger>
void create_clustering_graph(Graph& K,
                             const double N,
                             const double m,
                             const double c,
                             const double diameter,
                             const std::string degree_distr,
                             const bool calculate_c,
                             const bool calculate_diam,
                             const double tolerance,
                             RNGType& rng,
                             std::uniform_real_distribution<double>& distr,
                             const Logger& log)
{

    log->info("Assembling clustering component ... ");

    // Target values
    const double N_target = (diameter > 0)
        ?  std::round(N/num_vertices(K))
        :  N;
    const double m_target = (diameter > 0)
        ? std::round((m+1)/(K[0].state.mean_deg+1)-1)
        : m;

    log->info("Clustering component target size: N={}, m={}", N_target, m_target);

    // Get factors of N_target, m_target
    std::vector<double> N_factors;
    std::vector<double> m_factors;
    for (const auto& factor : Utils::closest_N_factors(N_target)) {
        N_factors.push_back(factor.first);
        N_factors.push_back(factor.second);
    }
    for (const auto& factor : Utils::closest_mean_deg_factors(m_target)) {
        m_factors.push_back(factor.first);
        m_factors.push_back(factor.second);
    }

    // ... Assembly start ......................................................

    // Component graphs & properties
    Graph G{}, H{};
    double N_G=N_target, N_H=1, m_G=m_target, m_H=0, var_G, var_H=0;
    double c_K, c_G, c_H=1, diam_G=-1, diam_H;

    // ... Base graph creation .................................................
    // Base graph is the graph generated using the target values for N
    if (degree_distr=="scale-free") {
        G = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N_target, m_target,
                                                              false, rng);
        c_G = global_clustering_coeff(G);
    }
    else {
        G = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N_target, m_target,
                                                          false, false, rng);
        c_G = m_target/(N_target-1); // this can cause problems ...
    }

    auto ds_G = degree_statistics(G);

    // If diameter is specified: c_G must include the clustering coefficient
    // from the previous assembly.
    if (diameter > 0) {
        c_G = Utils::Kronecker_clustering(c_G, K[0].state.clustering_global,
                                          ds_G.first, K[0].state.mean_deg,
                                          ds_G.second, K[0].state.var);
        var_G = Utils::Kronecker_degree_variance(ds_G.first, K[0].state.mean_deg,
                                                 ds_G.second, K[0].state.var);
        diam_G = std::max(diameter, Utopia::Models::NetworkAnalyser::diameter(G));
    }
    else {
        var_G = ds_G.second;
    }

    // Relative errors of base case
    double error = err_func(
        {rel_err(ds_G.first, m_target),
         rel_err(c_G, c),
         diameter_error(diam_G, diameter)}
    );

    // Output base graph property info
    log->info("Base graph created: N_G={}, m_G={}, c_G={}, {}error={}",
               N_G, m_G, c_G,
               (diameter > -1 ? "diam_G="+to_string(diam_G)+", " : ""),
               error);

    // ... Grid search .........................................................
    // If the base case is insufficent: grid search over N_fac, m_fac

    // Whether second component H has clustering coefficient 0
    bool zero_c = false;

    if (error > tolerance) {
        log->debug("Commencing grid search ...");
        // Temporary graph and properties
        Graph T{};
        double m_T, c_T, diam_T=-1, var_T;
        // Temporary properties of graph H
        double N_H_temp, m_H_temp, c_H_temp, diam_H_temp=-1;
        // Predicted c_K of resulting Kronecker product
        double predicted_c_K;
        // Current error value
        double current_err;

        for (const auto& m_fac: m_factors) {

            if (degree_distr=="scale-free" and m_fac < 3) { continue; }

            for (const auto& n_fac : N_factors) {

                if (m_fac >= n_fac) { continue; }

                log->debug("Checking case N_G={}, m_G={} ... ", n_fac, m_fac);

                // Generate temporary graph from current factorisation of N, m
                if (degree_distr=="scale-free") {
                    T = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(
                        n_fac, m_fac, false, rng);
                    c_T = global_clustering_coeff(T);
                }
                else {
                    T = Utopia::Graph::create_ErdosRenyi_graph<Graph>(
                        n_fac, m_fac, false, false, rng);
                    c_T = m_fac/(n_fac-1); // this can cause problems ...
                }

                const auto ds_T = degree_statistics(T);

                if (diameter > 0) {
                    m_T = Utils::Kronecker_mean_degree(ds_T.first, K[0].state.mean_deg);
                    var_T = Utils::Kronecker_degree_variance(ds_T.first,
                                                             K[0].state.mean_deg,
                                                             ds_T.second,
                                                             K[0].state.var);
                    diam_T = std::max(diameter, Utopia::Models::NetworkAnalyser::diameter(T));

                    c_T = Utils::Kronecker_clustering(c_T, K[0].state.clustering_global,
                                                ds_T.first, K[0].state.mean_deg,
                                                ds_T.second, K[0].state.var);
                }
                else {
                    m_T = ds_T.first;
                    var_T = ds_T.second;
                }

                // Find suitable component H
                for (c_H_temp = 0; c_H_temp < 2; ++c_H_temp) {

                    for (int a = 0; a < 2; ++a) {
                        if (c_H_temp == 0 and a == 1) {continue;}
                        // First attempt: use calculated m_H_temp from get_mean_deg_c
                        if (a==0) {
                            m_H_temp = std::round(std::max(
                                  2.,
                                  Utils::get_mean_deg_c(c_T, c, m_T, var_T)));
                            N_H_temp = (c_H_temp == 1)
                                ? m_H_temp+1
                                : std::max(m_H_temp+1, std::round(N_target/n_fac));
                        }
                        // Second attempt: Use the complementary value to m_fac
                        // instead of the calculated mean degree
                        else {
                            m_H_temp = std::round((m_target+1)/(m_fac+1))-1;
                            N_H_temp = N_target/n_fac;
                        }

                        if (m_H_temp >= N_H_temp) { continue; }

                        // Calculate predicted clustering coefficient (var_H = 0)
                        predicted_c_K = Utils::Kronecker_clustering(
                              c_T, c_H_temp, m_T, m_H_temp, var_T, var_H);

                        // Estimate diameter of component H
                        if (diameter > 0) {
                            diam_H_temp = N_H_temp/m_H_temp;
                        }

                        // Calculate error
                        current_err = err_func(
                            {rel_err(n_fac*N_H_temp, N_target),
                             rel_err((ds_T.first+1)*(m_H_temp+1), m_target+1),
                             rel_err(predicted_c_K, c),
                             diameter_error(std::max(diam_T, diam_H_temp), diameter)});

                        // If the current graph reduces the error: set graph factor
                        // values to current state in order to later reproduce the
                        // graph T; the two generation values (N_G and m_G) do not
                        // include the diameter component, though the topology
                        // values (diamt_T, c_T and var_T) do
                        if (current_err < error) {
                            // Generation values for G
                            N_G=n_fac;
                            m_G=std::round(ds_T.first);
                            // Generation values for H
                            N_H=N_H_temp;
                            m_H=m_H_temp;
                            // Topology values for G
                            c_G=c_T;
                            diam_G=diam_T;
                            var_G=ds_T.second;
                            // Topology values for H
                            c_H = c_H_temp;
                            // Predicted c_K
                            c_K=predicted_c_K;
                            // Whether or not H has clustering 0
                            zero_c = !(c_H==1 && a==0);
                            // Set the error to the current value
                            error = current_err;
                            log->debug("Improvement: N_T={}, m_T={}, c_T={}, "
                              "diam_T={}, N_H={}, m_H={}, c_H={}, predicted c_K={}, "
                              "error={}. {}",
                              N_G, m_G, c_G, diam_G, N_H, m_H, c_H,
                              predicted_c_K, error);
                        }
                    }

                }
                if (error < tolerance) {break;}
            }
            if (error < tolerance) {break;}
        }

        log->debug("Grid search complete.");
    }

    // If base graph is sufficient: no need to regenerate the base graph
    if (N_H == 1) {
        log->info("Base graph G is sufficient");
    }

    else {
        log->info("Generating graph G ...");

        // Generate the components G and H based on the values obtained from the
        // grid search
        if (degree_distr == "scale-free") {
            G = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N_G, m_G,
                                                                  false, rng);
        }
        else {
            G = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N_G, m_G,
                                                              false, false, rng);
        }
        log->info("Done. Generating graph H ...");
    }

    Utils::add_self_edges(G);


    if (N_H == 1) {
        H = Graph{1};
    }

    // H is a zero-clustering graph
    else if (zero_c) {
        if (N_H >= 2.*m_H) {
            // To do: Adapt this properly to m
            if (static_cast<int>(N_H) % 2) {
                ++N_H;
            }
            H = AuxGraphs::create_zero_c_graph<Graph>(N_H, m_H);
        }
        else {
            if (static_cast<int>(m_H) % 2) {
              ++m_H;
            }
            H = Utopia::Graph::create_regular_graph<Graph>(N_H, m_H, false);
        }
    }

    // H is a complete graph
    else {
        H = Utopia::Graph::create_complete_graph<Graph>(m_H+1);
    }

    Utils::add_self_edges(H);

    // Calculate properties
    if (diameter > 0 or calculate_diam) {
        diam_G = std::max(diameter, Utopia::Models::NetworkAnalyser::diameter(G));
        diam_H = Utopia::Models::NetworkAnalyser::diameter(H);
    }

    log->info("Done: Results: N_G={}, m_G={}, c_G={}, N_H={}, m_H={}, c_H={};"
              " predicted c_K={}{}",
              N_G, m_G, c_G, N_H, m_H, c_H, c_K,
              (diameter > 1)
                ? "; diam_G="+to_string(diam_G)+", diam_H="+to_string(diam_H)
                : "");

    // Combine G with graph from previous assembly, if given
    if (diameter > 1){
        log->info("Kronecker product of G with component from previous assembly ...");
        G = Utils::Kronecker_product(K, G, rng, distr);
    }

    // ... Create Kronecker graph and write properties .........................
    if (N_H == 1) {
        K = G;
    }
    else {
        K = Utils::Kronecker_product(G, H, rng, distr);
    }

    if (calculate_c) {
        if (num_vertices(H)==1) {
            K[0].state.clustering_global = c_G;
        }
        else {
            K[0].state.clustering_global =
                Utils::Kronecker_clustering(c_G, c_H, m_G, m_H, var_G, var_H);
        }
        log->info("Kronecker product clustering_coeff: {}", K[0].state.clustering_global);
    }

    if (calculate_diam) {
        K[0].state.diameter = std::max(diam_G, diam_H);
        log->info("Kronecker product diameter: {}", K[0].state.diameter);
    }
}

}
#endif