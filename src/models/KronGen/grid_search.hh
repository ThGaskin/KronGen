#ifndef UTOPIA_MODELS_KRONGEN_GRIDSEARCH
#define UTOPIA_MODELS_KRONGEN_GRIDSEARCH

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

#include "type_definitions.hh"
#include "objective_functions.hh"
#include "utils.hh"

namespace Utopia::Models::KronGen::GridSearch {

using namespace Utopia::Models::KronGen;
using namespace Utopia::Models::KronGen::TypeDefinitions;

// ... Grid search utility functions ...........................................

// Check function for invalid factor combinations: number of vertices must be
// larger than the mean degree, and mean degree 1 is only permissible if N = 2
template<typename TargetMap>
bool check_validity(const factor& n,
                    const factor& k,
                    TargetMap& analysis_and_targets)
{
    bool is_scale_free = Utils::is_param(analysis_and_targets, "degree_sequence");

    if (n.size() != k.size()) {
        return false;
    }
    if (n.size() == 0 or k.size() == 0) {
        return false;
    }

    for (size_t i = 0; i < k.size(); ++i) {

        // mean degrees must be smaller than number of vertices
        if (k[i]>=n[i]){
            return false;
        }

        if (k[i]<=1 and n[i]!=(k[i]+1)) {
            return false;
        }
        if (is_scale_free) {
            if (k[i] < 3 and n[i] != (k[i]+1)) {
                return false;
            }
        }
    }

    return true;
}

// Find positions for possible graph types
template<typename TargetMap>
std::vector<std::vector<GraphType>> find_possible_types (
                          const factor& n,
                          const factor& k,
                          TargetMap& analysis_and_targets,
                          const double& grid_size)
{

    std::vector<std::vector<GraphType>> type_list
    = Utils::is_param(analysis_and_targets, "degree_sequence")
      ? std::vector<std::vector<GraphType>>{{n.size(), GraphType::KlemmEguiluz}}
      : std::vector<std::vector<GraphType>>{{n.size(), GraphType::ErdosRenyi}};

    // Include chain graphs if the diameter is a target parameter
    if (Utils::is_param(analysis_and_targets, "diameter")) {
        double err = grid_size;
        int best_candidate = -1;
        for (size_t i = 0; i < n.size(); ++i) {
            if (k[i] > 2) {
                continue;
            }
            double curr_est
            = ObjectiveFuncs::err_func(n[i]-1,
                std::any_cast<double>(analysis_and_targets["diameter"]["target"].second));

            if ((curr_est < err)) {
                best_candidate = i;
                err = curr_est;
            }
        }

        if (best_candidate != -1) {
            auto new_type = type_list[0];
            new_type[best_candidate] = GraphType::Chain;

            type_list.emplace_back(new_type);
        }

    }

    // Include regular graphs
    if (!Utils::is_param(analysis_and_targets, "degree_sequence")) {
        const auto g = type_list.size();
        for (size_t t = 0; t < g; ++t) {
            int best_candidate = -1;
            for (size_t i = 0; i< k.size(); ++i){
                if ((k[i]%2 == 0))
                {
                    best_candidate = i;
                    auto new_type = type_list[t];
                    new_type[i] = GraphType::Regular;
                    type_list.emplace_back(new_type);
                }
            }
            if (best_candidate != -1) {
                auto new_type = type_list[t];
                new_type[best_candidate] = GraphType::Regular;
                type_list.emplace_back(new_type);
            }
        }
    }

    return type_list;
}

// Generates a grid of graph descriptors
template <typename TargetMap, typename Logger, typename RNGType>
ParetoSet generate_grid(const std::pair<size_t, size_t>& grid_center,
                        const double& grid_size,
                        const size_t& min_factors,
                        const size_t& max_factors,
                        TargetMap& analysis_and_targets,
                        const Logger& log,
                        RNGType& rng)
{

    ParetoSet grid = {};

    std::uniform_real_distribution<double> prob_distr;

    log->debug("Factorising N and k ...");

    auto n_grid = Utils::get_N_grid(grid_center.first, grid_size,
                                    min_factors, max_factors);
    auto k_grid = Utils::get_k_grid(grid_center.second, grid_size,
                                    min_factors, max_factors);

    // If no factorisation is possible, use grid centers
    if (n_grid.empty() or k_grid.empty()) {
        log->warn("No factorisation of either N or k possible! Increase values.");
        n_grid = {{grid_center.first}};
        k_grid = {{grid_center.second}};
    }

    double mu_0 = prob_distr(rng);
    if (Utils::is_param(analysis_and_targets, "clustering")) {
        double c_T = std::any_cast<double>(analysis_and_targets["clustering"]["target"].second);
        mu_0 = Utils::mu_inv(grid_center.first, grid_center.second, c_T);
    }

    for (const auto& n : n_grid) {

        for (const auto& k : k_grid) {

            // Ascertain n and k factors are valid
            if (not check_validity(n, k, analysis_and_targets)) {
                continue;
            }

            // Collect possible graph types
            const auto types = find_possible_types(n, k,
                analysis_and_targets, grid_size);

            for (const auto& t : types) {
                ParetoPoint p = {};

                for (size_t i = 0; i < n.size(); ++ i){
                    const GraphDesc g = {n[i], k[i], t[i], mu_0};
                    p.emplace_back(g);
                }
                grid.emplace_back(p);
            }
        }
    }

    return grid;
}

// ... Grid search functions ...................................................

/// AOF grid search: return candidates from a grid that minimise objective an
/// aggregate Hamiltonian.
/**
  * \tparam Logger              The logger type
  *
  * \param min_factors          Minimum number of Kronecker factors
  * \param max_factors          Maximum number of Kronecker factors
  * \param grid_center          Center of the grid to be searched
  * \param grid_size            The size of the grid
  * \param targets              The names and values of the target parameters
  * \param objective_funcs      The weights and objective functions used in the
                                aggregate objective function
  * \param log                  The model logger
  *
  * \return ParetoSet           A vector of graph tuples (Pareto points)
  */
template<typename Logger, typename RNGType, typename TargetMap, typename ObjFuncMap>
ParetoSet grid_search(const size_t& min_factors,
                            const size_t& max_factors,
                            const std::pair<size_t, size_t>& grid_center,
                            const double& grid_size,
                            TargetMap& analysis_and_targets,
                            ObjFuncMap& objective_funcs,
                            const Logger& log,
                            RNGType& rng)
{

    // Get the graph grid
    ParetoSet Grid = generate_grid(
      grid_center, grid_size, min_factors, max_factors,
      analysis_and_targets, log, rng
    );

    // Initialise vector for resulting factors and resulting error
    ParetoSet Paretos = {};
    double global_error = 0;

    // Extract target parametes from the analysis_and_targets map and place into
    // separate map
    const auto targets = Utils::extract_targets(analysis_and_targets);

    // Perform grid search
    log->debug("Commencing grid search ...");

    for (const auto& graphs : Grid) {

        // Calculate error function for the list of graphs considered
        double current_err = 0;
        for (const auto& p : targets) {

            const std::string param = p.first;

            const double weight = objective_funcs[param].first;
            const auto err = objective_funcs[param].second(p.second, graphs, rng);
            current_err += weight*err;
        }

        // If error is equal or reduced, add current factorisation
        // to set of Pareto points. If error is reduced, delete
        // previous Pareto points
        if ((current_err <= global_error) or (Paretos.empty())){
            if (current_err < global_error) {
                Paretos = {};
            }
            Paretos.emplace_back(graphs);
            global_error = current_err;
        }
    }

    // Write the global error value and return the Pareto points
    analysis_and_targets["error"]["calculate"].second = global_error;

    return Paretos;
}

// Grid search for scale-free graphs (WIP)
template<typename Graph, typename Logger, typename RNGType, typename TargetMap>
ParetoSet grid_search_SF(const size_t& num_factors,
                               const std::pair<size_t, size_t>& grid_center,
                               const double grid_size,
                               TargetMap& analysis_and_targets,
                               const Logger& log,
                               RNGType& rng)
{
    // TO DO: This must be generalised to arbitrary target parameters
    const double c_T = std::any_cast<double>(analysis_and_targets["clustering"]["target"].second);

    // Get the graph grid
    // TO DO: Setting grid size in N, k to 0 returns empty or single grid
    // when N or k have no factors. Example: k = 150
    ParetoSet Grid = generate_grid(
      grid_center, 0, num_factors, num_factors, analysis_and_targets, log, rng
    );

    // Initialise vector for resulting factors and resulting error
    ParetoSet Paretos = {};

    // Initialise optimisation error and get target parameter
    double global_error = 0;

    // Get estimated ideal mu parameter from grid
    double mu_0 = Grid[0][0].mu;

    // Perform grid search
    log->debug("Commencing grid search ...");
    for (const auto& g : Grid) {

        auto graphs = g;
        // Take out largest factor in graph list
        size_t n_max = 0, n_max_idx = 0;
        for (size_t i = 0; i < graphs.size(); ++i) {
            if (graphs[i].num_vertices > n_max){
                n_max = graphs[i].num_vertices;
                n_max_idx = i;
            }
        }
        graphs.erase(graphs.begin(), graphs.begin()+n_max_idx);

        // Generate Kronecker graph from all factors except the largest component
        // and calculate error
        double n_G = 1;
        double k_G = 0, c_G = 0, var_G = 0;

        // TO DO: Generalise this for arbitrary target parameters!
        for (const auto& graph : graphs) {

            const Graph H = AuxGraphs::create_graph<Graph>(graph, rng);
            const double c_H = NetworkAnalyser::global_clustering_coeff(H);
            const auto deg_stats = NetworkAnalyser::degree_statistics(H);

            c_G = Utils::Kronecker_clustering(c_H, c_G, deg_stats.first, k_G,
                                              deg_stats.second, var_G);

            var_G = Utils::Kronecker_degree_variance(deg_stats.first, k_G,
                                                     deg_stats.second, var_G);

            k_G = Utils::Kronecker_mean_degree(deg_stats.first, k_G);
            n_G = Utils::Kronecker_num_vertices(graph.num_vertices, n_G);
        }

        // TO DO: control the fineness of mu grid from frontend
        const auto mu_grid = Utils::get_mu_grid(mu_0, grid_size, 10);

        // Perform a grid search over mu for the largest component
        auto special_graph = g[n_max_idx];
        for (const auto& mu : mu_grid) {
            special_graph.mu = mu;

            // Generate graph
            const Graph H = AuxGraphs::create_graph<Graph>(special_graph, rng);
            const double c_H = NetworkAnalyser::global_clustering_coeff(H);
            const auto deg_stats = NetworkAnalyser::degree_statistics(H);

            // Calculate resulting properties
            double n_res =  Utils::Kronecker_num_vertices(special_graph.num_vertices, n_G);
            double k_res =  Utils::Kronecker_mean_degree(deg_stats.first, k_G);
            double c_res =  Utils::Kronecker_clustering(c_H, c_G,
                                                        deg_stats.first, k_G,
                                                        deg_stats.second, var_G);

            // Calculate resulting error
            double current_error = (ObjectiveFuncs::err_func(n_res, grid_center.first)
                                + ObjectiveFuncs::err_func(k_res, grid_center.second)
                                + ObjectiveFuncs::err_func(c_res, c_T));

            // If error is equal or reduced, add current factorisation
            // to set of Pareto points. If error is reduced, delete
            // previous Pareto points
            if ((current_error <= global_error) or (Paretos.empty())){

                if (current_error  < global_error) {
                    Paretos = {};
                }

                Paretos.emplace_back(graphs);
                Paretos.back().emplace_back(special_graph);

                global_error = current_error;
            }
        }
        break;
    }

    // Write the global error value and return the Pareto points
    analysis_and_targets["error"]["calculate"].second = global_error;

    return Paretos;
}

// .............................................................................
// Get Pareto points using appropriate grid search function
template<typename Graph, typename Logger, typename RNGType, typename TargetMap, typename ObjFuncMap>
ParetoSet get_Paretos(const size_t& min_factors,
                      const size_t& max_factors,
                      const std::pair<size_t, size_t>& grid_center,
                      const double& grid_size,
                      TargetMap& analysis_and_targets,
                      ObjFuncMap& objective_funcs,
                      const Logger& log,
                      RNGType& rng)
{
    if (Utils::is_target(analysis_and_targets, "degree_sequence", "scale-free")) {
        return grid_search_SF<Graph>(min_factors, grid_center, grid_size,
                                     analysis_and_targets, log, rng);
    }
    else {
        return grid_search(min_factors, max_factors, grid_center, grid_size,
                           analysis_and_targets, objective_funcs, log, rng);
    }
}

} // namespace KronGen::GridSearch

#endif // UTOPIA_MODELS_KRONGEN_GRIDSEARCH
