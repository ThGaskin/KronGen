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
    Graph G{}, H{};
    double N_G, N_H=1, m_G, m_H=0, var_G, var_H, c_K, c_G, c_H, diam_G, diam_H;
    const double tol = 0.05;

    if (diameter > 0) {
        double N_GG, m_GG, c_GG;
        N_GG = std::round(N/num_vertices(K));
        m_GG = std::round((m+1)/(K[0].state.mean_deg+1)-1);

        log->info("Clustering component target size: N={}, m={}", N_GG, m_GG);

        // Get factors of N, m .................................................
        std::vector<double> N_fac;
        std::vector<double> m_fac;
        for (const auto& n : Utils::closest_N_factors(N_GG)) {
            N_fac.push_back(n.first);
            N_fac.push_back(n.second);
        }
        for (const auto& mm : Utils::closest_mean_deg_factors(m_GG)) {
            m_fac.push_back(mm.first);
            m_fac.push_back(mm.second);
        }
        log->info("Have {} possible factors for N, {} possible factors for m.", N_fac.size(), m_fac.size());

        // Base case: one Kronecker product .....................................

        Graph GG{};
        if (degree_distr=="scale-free") {
            GG = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N_GG, m_GG, false, rng);
            c_GG = global_clustering_coeff(GG);
        }
        else {
            GG = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N_GG, m_GG, false, false, rng);
            c_GG = m/(N_GG-1); // can cause issues
        }

        N_G = N;
        m_G = m;
        auto ds_GG = degree_statistics(GG);
        c_G = Utils::Kronecker_clustering(c_GG, K[0].state.clustering_global,
                                   ds_GG.first, K[0].state.mean_deg,
                                   ds_GG.second, K[0].state.var);

        var_G = Utils::Kronecker_degree_variance(ds_GG.first, K[0].state.mean_deg,
                                          ds_GG.second, K[0].state.var);
        diam_G = std::max(diameter, Utopia::Models::NetworkAnalyser::diameter(GG));

        double diam_err = (diam_G > diameter)
            ? Utils::rel_err(diam_G, diameter)
            : 0;
        double error = Utils::err_func(
            {Utils::rel_err(ds_GG.first, m), Utils::rel_err(c_G, c), diam_err}
        );
        log->info("Base graph created: N_G={}, m_G={}, c_G={}, diam_G={}, error={}",
                   N_G, m_G, c_G, diam_G, error);

        // ... Grid search .....................................................
        bool zero_c;

        // If the base case is insufficient: grid search over N_fac, m_fac .....
        if (error > tol) {
            for (const auto& mm_G: m_fac) {

                if (degree_distr=="scale-free" and mm_G < 3) { continue; }

                for (const auto& nn_G : N_fac) {

                    if (mm_G >= nn_G) { continue; }

                    log->info("Checking N_G={}, m_G={} ... ", nn_G, mm_G);

                    // Temporary parameters
                    double cc_K, cc_H=0, ddiam_G, nn_H, mm_H, vv_H=0, err;

                    // Generate temporary graph
                    if (degree_distr=="scale-free") {
                        GG = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(
                            nn_G, mm_G, false, rng
                        );
                        c_GG = global_clustering_coeff(GG);
                    }
                    else {
                        GG = Utopia::Graph::create_ErdosRenyi_graph<Graph>(
                            nn_G, mm_G, false, false, rng
                        );
                        c_GG = (mm_G)/(nn_G-1); // this can cause problems ...
                    }
                    ds_GG = degree_statistics(GG);
                    ddiam_G = std::max(diameter, Utopia::Models::NetworkAnalyser::diameter(GG));
                    c_GG = Utils::Kronecker_clustering(c_GG, K[0].state.clustering_global,
                                                ds_GG.first, K[0].state.mean_deg,
                                                ds_GG.second, K[0].state.var);

                    // ... First attempt: Calculated mean degree ...............
                    mm_H = std::round(std::max(2.,
                      Utils::get_mean_deg_c(c_GG, c, ds_GG.first, ds_GG.second)));

                    for (cc_H=0; cc_H<2; ++cc_H) {

                        // Get number of vertices
                        nn_H = (cc_H == 1) ? mm_H+1 : std::max(mm_H+1, N_GG/nn_G);

                        // Predicted clustering coefficient
                        cc_K = Utils::Kronecker_clustering(
                              c_GG, cc_H, ds_GG.first, mm_H, ds_GG.second, vv_H
                        );

                        diam_err = (std::max(ddiam_G, nn_H/mm_H) > diameter)
                            ? Utils::rel_err(std::max(ddiam_G, nn_H/mm_H), diameter)
                            : 0;

                        // Calculate error
                        err = Utils::err_func({
                            Utils::rel_err(nn_G*nn_H, N_GG),
                            Utils::rel_err((ds_GG.first+1)*(mm_H+1), m_GG+1),
                            Utils::rel_err(cc_K, c),
                            diam_err
                        });

                        // If current graph minimises error: set graph factor
                        // values to current state
                        if (err < error) {
                            N_G=nn_G;
                            N_H=std::round(nn_H);
                            m_G=mm_G;
                            m_H=mm_H;
                            c_K=cc_K;
                            c_G=c_GG;
                            diam_G = std::max(ddiam_G, nn_H/mm_H);
                            var_G=ds_GG.second;
                            zero_c = (cc_H == 0);
                            error = err;
                        }
                    }

                    // ... Second attempt: inverse mean degree .................
                    mm_H = std::round((m_GG+1)/(mm_G+1))-1;
                    nn_H = N_GG/nn_G;
                    if (mm_H >= nn_H) {continue;}
                    for (cc_H=0; cc_H<2; ++cc_H) {
                        cc_K = Utils::Kronecker_clustering(c_GG, cc_H, ds_GG.first, mm_H, ds_GG.second, vv_H);

                        diam_err = (std::max(ddiam_G, nn_H/mm_H) > diameter)
                            ? Utils::rel_err(std::max(ddiam_G, nn_H/mm_H), diameter)
                            : 0;

                        err = Utils::err_func(
                        {
                            Utils::rel_err(nn_G*nn_H, N_GG),
                            Utils::rel_err((ds_GG.first+1)*(mm_H+1), m+1),
                            Utils::rel_err(cc_K, c),
                            diam_err
                        });

                        if (err < error) {
                            N_G=nn_G;
                            N_H=std::round(nn_H);
                            m_G=mm_G;
                            m_H=mm_H;
                            c_K=cc_K;
                            c_G=c_GG;
                            var_G=ds_GG.second;
                            zero_c = (cc_H == 0);
                            diam_G = std::max(ddiam_G, nn_H/mm_H);
                            error = err;
                        }
                    }
                    if (error < tol) {break;}
                }
                if (error < tol) {break;}
            }
        }
        // ... End of grid search ..............................................
        log->info("Grid search complete. Generating graphs G and H ...");
        if (degree_distr == "scale-free") {
            G = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N_G, m_G, false, rng);
        }
        else {
            G = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N_G, m_G, false, false, rng);
        }

        if (N_G == N or m_G == m) {
            H = Utopia::Graph::create_complete_graph<Graph>(1);
        }

        else if (zero_c) {
            if (N_H >= 2.*m_H) {
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
        else {
            H = Utopia::Graph::create_complete_graph<Graph>(N_H);
        }

        log->info("Done: Results: N_G={}, m_G={}, N_H={}, m_H={}, c_K={}", num_vertices(G), m_G, N_H, m_H, c_K);

        if (calculate_diam) {
            const auto s = fourSweep<vertices_size_type>(G);
            diam_G = iFUB(s.first, s.second, 0, G);
        }

        Utils::add_self_edges(G);
        G = Utils::Kronecker_product(K, G, rng, distr);

        log->info("Kronecker product with diameter component complete; N_G={}; "
                  "factor properties: c_G={}, diam_G={}",
                   num_vertices(G), c_G, diam_G);

    }

                                                                                else {

                                                                                    // Get factors of N, m .................................................
                                                                                    std::vector<double> N_fac;
                                                                                    std::vector<double> m_fac;
                                                                                    const auto N_temp = Utils::closest_N_factors(N);
                                                                                    for (const auto& n : N_temp) {
                                                                                        N_fac.push_back(n.first);
                                                                                        N_fac.push_back(n.second);
                                                                                    }
                                                                                    const auto m_temp = Utils::closest_mean_deg_factors(m);
                                                                                    for (const auto& mm : m_temp) {
                                                                                        m_fac.push_back(mm.first);
                                                                                        m_fac.push_back(mm.second);
                                                                                    }

                                                                                    // Base case: no Kronecker product .....................................
                                                                                    Graph GG{};
                                                                                    if (degree_distr=="scale-free") {
                                                                                        GG = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N, m, false, rng);
                                                                                        c_G = global_clustering_coeff(GG);
                                                                                    }
                                                                                    else {
                                                                                        GG = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N, m, false, false, rng);
                                                                                        c_G = m/(N-1); // this can cause problems ...
                                                                                    }
                                                                                    auto ds_G = degree_statistics(GG);
                                                                                    N_G = N;
                                                                                    m_G = m;
                                                                                    var_G = ds_G.second;
                                                                                    double error = Utils::err_func(
                                                                                        {Utils::rel_err(ds_G.first, m), Utils::rel_err(c_G, c)}
                                                                                    );
                                                                                    bool zero_c;

                                                                                    // ... Grid search .....................................................
                                                                                    // If the base case is insufficient: grid search over N_fac, m_fac .....
                                                                                    if (error > tol) {
                                                                                        for (const auto& mm_G: m_fac) {

                                                                                            if (degree_distr=="scale-free" and mm_G < 3) { continue; }

                                                                                            for (const auto& nn_G : N_fac) {

                                                                                                if (mm_G >= nn_G) { continue; }

                                                                                                log->info("Checking N_G={}, m_G={} ... ", nn_G, mm_G);

                                                                                                // Temporary parameters
                                                                                                double cc_G, cc_K, cc_H=0, nn_H, mm_H, vv_H=0, err;

                                                                                                // Generate temporary graph
                                                                                                if (degree_distr=="scale-free") {
                                                                                                    GG = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(
                                                                                                        nn_G, mm_G, false, rng
                                                                                                    );
                                                                                                    cc_G = global_clustering_coeff(GG);
                                                                                                }
                                                                                                else {
                                                                                                    GG = Utopia::Graph::create_ErdosRenyi_graph<Graph>(
                                                                                                        nn_G, mm_G, false, false, rng
                                                                                                    );
                                                                                                    cc_G = mm_G/(nn_G-1);
                                                                                                }
                                                                                                ds_G = degree_statistics(GG);

                                                                                                // ... First attempt: Calculated mean degree ...............
                                                                                                mm_H = std::round(std::max(2.,
                                                                                                  Utils::get_mean_deg_c(cc_G, c, ds_G.first, ds_G.second)));

                                                                                                for (cc_H=0; cc_H<2; ++cc_H) {

                                                                                                    // Get number of vertices
                                                                                                    nn_H = (cc_H == 1) ? mm_H+1 : std::max(mm_H+1, N/nn_G);

                                                                                                    // Predicted clustering coefficient
                                                                                                    cc_K = Utils::Kronecker_clustering(
                                                                                                          cc_G, cc_H, ds_G.first, mm_H, ds_G.second, vv_H
                                                                                                    );

                                                                                                    // Calculate error
                                                                                                    err = Utils::err_func({
                                                                                                        Utils::rel_err(nn_G*nn_H, N),
                                                                                                        Utils::rel_err((ds_G.first+1)*(mm_H+1), m+1),
                                                                                                        Utils::rel_err(cc_K, c)
                                                                                                    });

                                                                                                    // If current graph minimises error: set graph factor
                                                                                                    // values to current state
                                                                                                    if (err < error) {
                                                                                                        N_G=nn_G;
                                                                                                        N_H=nn_H;
                                                                                                        m_G=mm_G;
                                                                                                        m_H=mm_H;
                                                                                                        c_K=cc_K;
                                                                                                        c_G=cc_G;
                                                                                                        var_G=ds_G.second;
                                                                                                        zero_c = (cc_H == 0);
                                                                                                        error = err;
                                                                                                    }
                                                                                                }

                                                                                                // ... Second attempt: inverse mean degree .................
                                                                                                mm_H = std::round((m+1)/(mm_G+1))-1;
                                                                                                nn_H = N/nn_G;
                                                                                                if (mm_H >= nn_H) {continue;}
                                                                                                for (cc_H=0; cc_H<2; ++cc_H) {
                                                                                                    cc_K = Utils::Kronecker_clustering(cc_G, cc_H, ds_G.first, mm_H, ds_G.second, vv_H);

                                                                                                    err = Utils::err_func(
                                                                                                    {
                                                                                                        Utils::rel_err(nn_G*nn_H, N),
                                                                                                        Utils::rel_err((ds_G.first+1)*(mm_H+1), m+1),
                                                                                                        Utils::rel_err(cc_K, c)
                                                                                                    });

                                                                                                    if (err < error) {
                                                                                                        N_G=nn_G;
                                                                                                        N_H=nn_H;
                                                                                                        m_G=mm_G;
                                                                                                        m_H=mm_H;
                                                                                                        c_K=cc_K;
                                                                                                        c_G=cc_G;
                                                                                                        var_G=ds_G.second;
                                                                                                        zero_c = (cc_H == 0);
                                                                                                        error = err;
                                                                                                    }
                                                                                                }
                                                                                                if (error < tol) {break;}
                                                                                            }
                                                                                            if (error < tol) {break;}
                                                                                        }
                                                                                    }
                                                                                    // ... End of grid search ..............................................
                                                                                    log->info("Done: Results: N_G={}, m_G={}, N_H={}, m_H={}, c_K={}", N_G, m_G, N_H, m_H, c_K);

                                                                                    if (degree_distr == "scale-free") {
                                                                                        G = Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N_G, m_G, false, rng);
                                                                                    }
                                                                                    else {
                                                                                        G = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N_G, m_G, false, false, rng);
                                                                                    }

                                                                                    if (N_G == N or m_G == m) {
                                                                                        H = Utopia::Graph::create_complete_graph<Graph>(1);
                                                                                    }

                                                                                    else if (zero_c) {
                                                                                        if (N_H >= 2.*m_H) {
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
                                                                                    else {
                                                                                        H = Utopia::Graph::create_complete_graph<Graph>(N_H);
                                                                                    }

                                                                                    if (calculate_diam) {
                                                                                        const auto s = fourSweep<vertices_size_type>(G);
                                                                                        diam_G = iFUB(s.first, s.second, 0, G);
                                                                                    }

                                                                                    log->info("First factor properties: c_G={}{}.", c_G,
                                                                                              (calculate_diam ? ", diam_G="+to_string(diam_G) : ""));

                                                                                    Utils::add_self_edges(G);

                                                                                }

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

    log->info("Kronecker product G x H; H properties: {}", num_vertices(H));

    K = Utils::Kronecker_product(G, H, rng, distr);

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
