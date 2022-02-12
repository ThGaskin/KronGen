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
    const bool is_scale_free = Utils::is_target(
        analysis_and_targets, "degree_sequence", "scale-free"
    );

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
                          const double& grid_range)
{

    std::vector<std::vector<GraphType>> type_list
    = Utils::is_param(analysis_and_targets, "degree_sequence")
      ? std::vector<std::vector<GraphType>>{{n.size(), GraphType::KlemmEguiluz}}
      : std::vector<std::vector<GraphType>>{{n.size(), GraphType::ErdosRenyi}};

    // Include chain graphs if the diameter is a target parameter
    if (Utils::is_param(analysis_and_targets, "diameter")) {
        double err = grid_range;
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
                        const double& grid_range,
                        const size_t& min_factors,
                        const size_t& max_factors,
                        TargetMap& analysis_and_targets,
                        const Logger& log,
                        RNGType& rng)
{

    ParetoSet grid = {};

    std::uniform_real_distribution<double> prob_distr;

    log->debug("Factorising N and k ...");

    auto n_grid = Utils::get_N_grid(grid_center.first, grid_range,
                                    min_factors, max_factors);
    auto k_grid = Utils::get_k_grid(grid_center.second, grid_range,
                                    min_factors, max_factors);

    // If no factorisation is possible, decrease grid center values
    if (n_grid.empty()) {
        log->warn("No factorisation of N = {} possible! Decreasing value.",
                  grid_center.first);
        const size_t min_n = pow(2, min_factors);
        std::pair<size_t, size_t> grid_center_mod = {std::max(min_n, grid_center.first-1), grid_center.second};
        return generate_grid(grid_center_mod, grid_range, min_factors, max_factors, analysis_and_targets, log, rng);
    }
    else if (k_grid.empty()) {
        log->warn("No factorisation of k = {} possible! Decreasing value.",
                  grid_center.second);
        std::pair<size_t, size_t> grid_center_mod = {grid_center.first, std::max(grid_center.second-1, size_t(5))};
        return generate_grid(grid_center_mod, grid_range, min_factors, max_factors, analysis_and_targets, log, rng);
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
                analysis_and_targets, grid_range);

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

    if (grid.empty()) {
        log->warn("Grid empty! Increase tolerance or change target values!");
        return {{GraphDesc{2, 1, GraphType::Complete, 0}}};
    }
    return grid;
}

// ... Grid search functions ...................................................

/// AOF grid search: return candidates from a grid that minimise objective an
/// aggregate Hamiltonian.
/**
  * \param min_factors          Minimum number of Kronecker factors
  * \param max_factors          Maximum number of Kronecker factors
  * \param grid_center          Center of the grid to be searched
  * \param grid_range           The range of the grid in N, k
  * \param targets              The names and values of the target parameters
  * \param objective_funcs      The weights and objective functions used in the
                                aggregate objective function
  * \param log                  The model logger
  * \param rng                  The model rng
  *
  * \return ParetoSet           A vector of graph tuples (Pareto points)
  */
template<typename Logger, typename RNGType, typename ObjFuncMap>
ParetoSet grid_search(const size_t& min_factors,
                      const size_t& max_factors,
                      const std::pair<size_t, size_t>& grid_center,
                      const double& grid_range,
                      map_type& analysis_and_targets,
                      ObjFuncMap& objective_funcs,
                      const Logger& log,
                      RNGType& rng)
{

    // Get the graph grid
    log->debug("Generating grid ...");
    ParetoSet Grid = generate_grid(
        grid_center, grid_range, min_factors, max_factors,
        analysis_and_targets, log, rng
    );

    // Initialise vector for optimal factors and global optimisation error
    ParetoSet Paretos = {};
    double global_error = 0;

    // Extract target parametes from the analysis_and_targets map and place into
    // separate map
    const auto targets = Utils::extract_targets(analysis_and_targets);

    // Perform grid search
    log->debug("Searching grid ...");
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
                Paretos.clear();
            }
            Paretos.emplace_back(graphs);
            global_error = current_err;
        }
    }

    // Write the global error value and return the Pareto points
    analysis_and_targets["error"]["calculate"].second = global_error;

    return Paretos;
}

// Grid search for scale-free graphs
/**
  * \param num_factors          Number of Kronecker factors
  * \param grid_center          Center of the grid to be searched
  * \param grid_range           The range of the grid for N, k, mu
  * \param num_mu               The size of mu_grid
  * \param analysis_and_targets The names and values of the target parameters
  * \param log                  The model logger
  * \param rng                  The model rng
  *
  * \return ParetoSet           A vector of graph tuples (Pareto points)
  */
template<typename Graph, typename Logger, typename RNGType>
ParetoSet grid_search_SF(const size_t& num_factors,
                         const std::pair<size_t, size_t>& grid_center,
                         const double grid_range,
                         const size_t& num_mu,
                         map_type& analysis_and_targets,
                         const Logger& log,
                         RNGType& rng)
{

    // Get the graph grid
    ParetoSet Grid = generate_grid(
        grid_center, 0, num_factors, num_factors, analysis_and_targets, log, rng
    );

    // Initialise vector for resulting factors and resulting error
    ParetoSet Paretos = {};
    double global_error = 0;

    // Get estimated ideal mu parameter from grid and generate mu grid
    const double mu_0 = Grid[0][0].mu;
    const auto mu_grid = Utils::get_mu_grid(mu_0, grid_range, num_mu);

    // Extract target parametes from the analysis_and_targets map and place into
    // separate map
    const auto targets = Utils::extract_targets(analysis_and_targets);

    // Perform grid search over N and k
    log->debug("Searching grid ...");
    for (const auto& g : Grid) {

        // Take out largest factor in graph list
        auto graphs = g;
        size_t n_max = 0, n_max_idx = 0;
        for (size_t i = 0; i < graphs.size(); ++i) {
            if ((1000 > graphs[i].num_vertices) and (graphs[i].num_vertices > n_max)){
                n_max = graphs[i].num_vertices;
                n_max_idx = i;
            }
        }
        graphs.erase(graphs.begin()+n_max_idx);

        // Create a map to keep track of analysis targets
        map_type analysis_values;
        for (const auto& p : targets) {
            if (p.first == "num_vertices") {
                analysis_values[p.first]
                = entry_type{{"calculate", std::make_pair<bool, double>(true, 1)}};
            }
            else {
                analysis_values[p.first]
                = entry_type{{"calculate", std::make_pair<bool, double>(true, 0)}};
            }
            if (p.first == "clustering") {
                analysis_values["degree_variance"]
                = entry_type{{"calculate", std::make_pair<bool, double>(true, 0)}};
                analysis_values["mean_degree"]
                = entry_type{{"calculate", std::make_pair<bool, double>(true, 0)}};
            }
        }

        // Generate Kronecker product of all graphs except factor selected for
        // grid search in mu
        for (const auto& graph : graphs) {

            const Graph H = AuxGraphs::create_graph<Graph>(graph, rng);

            Utils::calculate_properties(H, analysis_values, log);

        }

        // Perform a grid search over mu for the largest component
        auto special_graph = g[n_max_idx];

        log->debug("Searching grid in mu over graph with N={}, k={} ... ",
            special_graph.num_vertices, special_graph.mean_degree);

        for (const auto& mu : mu_grid) {
            special_graph.mu = mu;
            map_type temporary_analysis_values = analysis_values;

            // Generate graph and calculate resulting properties
            const Graph H = AuxGraphs::create_graph<Graph>(special_graph, rng);
            Utils::calculate_properties(H, temporary_analysis_values, log);

            // Calculate optimisation error. If error is equal or reduced,
            // add current factorisationto set of Pareto points. If error is
            // reduced, delete previous Pareto points
            double current_error = 0;
            for (const auto& p : targets) {
                current_error += ObjectiveFuncs::err_func(std::any_cast<double>(temporary_analysis_values[p.first]["calculate"].second),
                                                          p.second);
            }
            if ((current_error <= global_error) or (Paretos.empty())){

                if (current_error  < global_error) {
                    Paretos.clear();
                }

                Paretos.emplace_back(graphs);
                Paretos.back().emplace_back(special_graph);

                global_error = current_error;
            }
        }
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
                      const size_t& num_mu,
                      const std::pair<size_t, size_t>& grid_center,
                      const double& grid_range,
                      TargetMap& analysis_and_targets,
                      ObjFuncMap& objective_funcs,
                      const Logger& log,
                      RNGType& rng)
{
    if (Utils::is_target(analysis_and_targets, "degree_sequence", "scale-free")) {
        return grid_search_SF<Graph>(min_factors, grid_center, grid_range, num_mu,
                                     analysis_and_targets, log, rng);
    }
    else {
        return grid_search(min_factors, max_factors, grid_center, grid_range,
                           analysis_and_targets, objective_funcs, log, rng);
    }
}

} // namespace KronGen::GridSearch

#endif // UTOPIA_MODELS_KRONGEN_GRIDSEARCH
