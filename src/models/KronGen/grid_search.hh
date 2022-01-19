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

#include "graph_types.hh"
#include "objective_functions.hh"
#include "utils.hh"

namespace Utopia::Models::KronGen::GridSearch {

using namespace Utopia::Models::KronGen;
using namespace Utopia::Models::KronGen::GraphTypes;

using factor = typename std::vector<std::size_t>;
using GraphListType
= typename std::vector<std::tuple<factor, factor, std::vector<GraphType>>>;
using target_map_type
= typename std::map<std::string, std::map<std::string, std::pair<bool, std::any>>>;

// Check function for invalid factor combinations: number of vertices must be
// larger than the mean degree, and mean degree 1 is only permissible if N = 2
bool check_validity(const factor& n,
                    const factor& k,
                    target_map_type& analysis_and_targets)
{
    bool is_scale_free = false;
    if (analysis_and_targets.find("degree_sequence") != analysis_and_targets.end()){
        if (analysis_and_targets["degree_sequence"]["target"].first) {
            is_scale_free = true;
        }
    }

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
std::vector<std::vector<GraphType>> find_possible_types (
                          const factor& n,
                          const factor& k,
                          target_map_type& analysis_and_targets,
                          const double& grid_size)
{

    std::vector<std::vector<GraphType>> type_list
    = std::vector<std::vector<GraphType>>{{n.size(), GraphType::ErdosRenyi}};
    if (analysis_and_targets.find("degree_sequence") != analysis_and_targets.end()){
        if (analysis_and_targets["degree_sequence"]["target"].first == true) {
            type_list
            = std::vector<std::vector<GraphType>>{{n.size(), GraphType::KlemmEguiluz}};
        }
    }

    // Include chain graphs if the diameter is a target parameter
    bool diameter_target = false;
    if (analysis_and_targets.find("diameter") != analysis_and_targets.end()){
        if (analysis_and_targets["diameter"]["target"].first == true) {
            diameter_target = true;
        }
    }

    if (diameter_target == true) {
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
    const auto g = type_list.size();
    for (size_t t = 0; t < g; ++t) {
        // size_t n_max = 0;
        // int best_candidate = -1;
        for (size_t i = 0; i< k.size(); ++i){
            if ((k[i]%2 == 0))
            {
                // best_candidate = i;
                // n_max = n[i];
                auto new_type = type_list[t];
                new_type[i] = GraphType::Regular;
                type_list.emplace_back(new_type);
            }
        }
        // if (best_candidate != -1) {
        //     auto new_type = type_list[t];
        //     new_type[best_candidate] = GraphType::Regular;
        //     type_list.emplace_back(new_type);
        // }
    }

    return type_list;
}


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
  * \return GraphListType       A vector of tuples of vectors representing
                                the factor properties (N, k, type)
  */
template<typename Logger, typename RNGType, typename obj_func_map_type>
GraphListType grid_search(
    const size_t& min_factors,
    const size_t& max_factors,
    const std::pair<size_t, size_t>& grid_center,
    const double& grid_size,
    target_map_type& analysis_and_targets,
    obj_func_map_type& objective_funcs,
    const Logger& log,
    RNGType& rng)
{

    // Get grids in N and k using interval bounds
    log->debug("Setting up grid ...");
    auto n_grid = Utils::get_N_grid(grid_center.first, grid_size,
                                          min_factors, max_factors);
    auto k_grid = Utils::get_k_grid(grid_center.second, grid_size,
                                          min_factors, max_factors);

    // Extract target parametes from the analysis_and_targets map and place into
    // separate map
    const auto targets = Utils::extract_targets(analysis_and_targets);


    // If the grid is empty, return grid center
    if (n_grid.empty() or k_grid.empty()) {
        log->warn("No factorisation of either N or k possible! Increase values.");
        n_grid = {{grid_center.first}};
        k_grid = {{grid_center.second}};
    }

    // Initialise vector for resulting factors and resulting error
    std::vector<std::tuple<factor, factor, std::vector<GraphType>>> Paretos = {};
    double error = -1;

    // Perform grid search
    log->debug("Commencing grid search ...");
    for (const auto& n: n_grid) {
         for (const auto& k: k_grid) {

            // Ascertain n and k factors are valid
            if (not check_validity(n, k, analysis_and_targets)) {
                continue;
            }

            // Collect possible graph types
            const auto graph_types = find_possible_types(n, k,
                                               analysis_and_targets, grid_size);

            // Calculate error function for all types and all factors considered
            for (const auto& t : graph_types) {
                double current_err = 0;
                for (const auto& p : targets) {

                    const std::string param = p.first;

                    const double weight = objective_funcs[param].first;
                    const auto err = objective_funcs[param].second(p.second, n, k, t, rng);
                    current_err += weight*err;
                }

                // If error constant, add current factorisation to set of Pareto points
                if (current_err == error) {
                    Paretos.emplace_back(std::make_tuple(n, k, t));
                }

                // If error reduced, use curent factorisation
                else if((current_err < error) or (error == -1)) {
                    Paretos = {{n, k, t}};
                    error = current_err;
                }
            }
        }
    }

    // Write the global error value and return the Pareto points
    analysis_and_targets["error"]["calculate"].second = error;

    return Paretos;
}

} // namespace KronGen::GridSearch

#endif // UTOPIA_MODELS_KRONGEN_GRIDSEARCH
