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

// Check function for invalid factor combinations: number of vertices must be
// larger than the mean degree, and mean degree 1 is only permissible if N = 2
bool check_validity(const factor& n, const factor& k)
{
    if (n.size() != k.size()) {
        return false;
    }

    for (size_t i = 0; i < k.size(); ++i) {

        // mean degrees must be smaller than number of vertices
        if (k[i]>=n[i]){
            return false;
        }

        if (k[i]<=1 and n[i]!=k[i]+1) {
            return false;
        }
    }

    return true;
}

// Find positions for possible graph types (ER, regular, chain)
void find_possible_types (const factor& n,
                          const factor& k,
                          const std::vector<double>& target_parameters,
                          std::vector<std::vector<GraphType>>& graph_types,
                          const double& max_error)
{
    double chain_err = max_error;

    // Include chain graphs
    auto g = graph_types.size();
    if (target_parameters.size() > 2 and target_parameters.back() > 1) {
        for (size_t t = 0; t < g; ++t) {
            int current_candidate = -1;
            for (size_t i = 0; i < k.size(); ++i) {
                if (k[i] <= 2) {
                    double curr_est = ObjectiveFuncs::err_func(n[i]-1, target_parameters.back());
                    if ((curr_est < chain_err)) {
                        current_candidate = i;
                        chain_err = curr_est;
                    }
                }
            }
            if (current_candidate != -1) {
                auto new_type = graph_types[t];
                new_type[current_candidate] = GraphType::Chain;
                graph_types.emplace_back(new_type);
            }
        }
    }

    // Include regular graphs
    g = graph_types.size();
    for (size_t t = 0; t < g; ++t) {
        size_t n_max = 0;
        int current_candidate = -1;
        for (size_t i = 0; i< k.size(); ++i){
            if ((k[i]%2 == 0) and (n[i] > n_max))
            {
                current_candidate = i;
                n_max = n[i];
            }
        }
        if (current_candidate != -1) {
            auto new_type = graph_types[t];
            new_type[current_candidate] = GraphType::Regular;
            graph_types.emplace_back(new_type);
        }
    }
}


/// AOF grid search: loop over N, return candidates that minimise objective
/// function. No need to generate complete graph
/**
  * \tparam Logger              The logger type
  *
  * \param num_factors          Number of Kronecker factors
  * \param grid_center          Centre of N-k-grid
  * \param target_parameters    The target parameters
  * \param obj_func             The objective functions for each parameter
  * \param weights              The weights used in the aggregate objective
                                function
  * \param max_error            The maximum error used to determine grid
                                search bounds
  * \param log                  The model logger
  *
  * \return result              A vector of tuples of vectors representing
                                the factor properties (N, k, type)
  * @throws invalid_argument    If target parameters are given and number of
                                target parameters does not match that of
                                objective functions or weights
  */
template<typename Logger>
std::pair<std::vector<std::tuple<factor, factor, std::vector<GraphType>>>, double> grid_search(
    const size_t& min_factors,
    const size_t& max_factors,
    const std::pair<size_t, size_t>& grid_center,
    const std::vector<double>& target_parameters,
    std::vector<std::function<double(double, factor, factor, std::vector<GraphType>)>>& obj_func,
    const std::vector<double>& weights,
    const double& max_error,
    const Logger& log)
{

    // If no target parameters are given, simply return starting values
    if (target_parameters.empty()) {
        return {{{{grid_center.first}, {grid_center.second}, {GraphType::ErdosRenyi}}}, 0.};
    }

    // Number of target parameters, objective functions, and weights must be equal
    if ((target_parameters.size() != obj_func.size())
       or (obj_func.size() != weights.size())
       or (weights.size() != target_parameters.size())) {
           throw std::invalid_argument("Number of target parameters, objective "
                                       "functions, and weights do not match!");
    }

    // Get grids in N and k using interval bounds
    log->debug("Setting up grid ...");
    const auto n_grid = Utils::get_N_grid(grid_center.first,
                                          max_error, min_factors, max_factors);
    const auto k_grid = Utils::get_k_grid(grid_center.second,
                                          max_error, min_factors, max_factors);
    log->debug("Grid set up.");

    if (n_grid.empty() or k_grid.empty()) {
        log->warn("No factorisation of either N or k possible! Increase values.");
        return {{{{grid_center.first}, {grid_center.second}, {GraphType::ErdosRenyi}}}, 0.};
    }

    // Initialise vector for resulting factors and resulting error
    std::vector<std::tuple<factor, factor, std::vector<GraphType>>> res = {};
    double error = -1;

    // Perform grid search
    for (const auto& n: n_grid) {
         for (const auto& k: k_grid) {

            // Ascertain n and k factors are valid
            if (not check_validity(n, k)) {
                continue;
            }

            // List possible type combinations
            std::vector<std::vector<GraphType>> graph_types;
            graph_types.emplace_back(n.size(), GraphType::ErdosRenyi);
            find_possible_types(n, k, target_parameters, graph_types, max_error);

            // Calculate error function for all types and all factors considered
            for (const auto& t : graph_types) {
                double current_err = 0;
                for (size_t i = 0; i < target_parameters.size(); ++i) {
                    current_err += weights[i]*obj_func[i](target_parameters[i], n, k, t);
                }

                // If error constant, add current factorisation to set of Pareto points
                if (current_err == error) {
                    res.emplace_back(std::make_tuple(n, k, t));
                }

                // If error reduced, use curent factorisation
                else if((current_err < error) or (error == -1)) {
                    res = {{n, k, t}};
                    error = current_err;
                }
            }
        }
    }

    // Return Pareto points and the global error
    return {res, error};
}

} // namespace KronGen::GridSearch

#endif // UTOPIA_MODELS_KRONGEN_GRIDSEARCH
