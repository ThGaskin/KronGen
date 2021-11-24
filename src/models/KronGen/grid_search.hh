#ifndef UTOPIA_MODELS_KRONGEN_GRIDSEARCH
#define UTOPIA_MODELS_KRONGEN_GRIDSEARCH

#include <vector>

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
#include "graph_types.hh"
#include "utils.hh"

#include "../NetworkAnalyser/graph_metrics.hh"

namespace Utopia::Models::KronGen::GridSearch {

using namespace std;
using namespace boost;
using namespace Utopia::Models::KronGen;
using namespace Utopia::Models::KronGen::GraphTypes;
using namespace Utopia::Models::KronGen::Utils;
using namespace Utopia::Models::NetworkAnalyser;

using factor = typename std::vector<std::size_t>;

// Check for invalid factor combinations
bool check_validity(const factor n, const factor k){

  if (n.size() != k.size()) {
      return false;
  }

  for (int i = 0; i < k.size(); ++i) {

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

void find_possible_types (const factor& n,
                          const factor& k,
                          const std::vector<double>& target_parameters,
                          std::vector<std::vector<GraphType>>& graph_types,
                          const double& max_error)
{
    std::vector<GraphType> types(n.size(), GraphType::ErdosRenyi);
    graph_types.emplace_back(types);
    double chain_err = max_error;

    // Include chain graphs
    if (target_parameters.back() > 1) {
        size_t current_candidate = -1;
        for (int i = 0; i < k.size(); ++i) {
            if (k[i] <= 2) {
                double curr_est = Utils::err_func(n[i]-1, target_parameters.back());
                if ((curr_est < chain_err)) {
                    current_candidate = i;
                    chain_err = curr_est;
                }
            }
        }
        if (current_candidate != -1) {
            types[current_candidate] = GraphType::Chain;
            graph_types.emplace_back(types);
        }
    }

    // Include regular graphs
    if ((target_parameters.size() == 3 and target_parameters.back() <= 1)
        or (target_parameters.size() == 4)) {

        const auto g = graph_types.size();
        for (int t = 0; t < g; ++t) {
            size_t n_max = 0;
            size_t current_candidate = -1;
            for (int i = 0; i< k.size(); ++i){
                if ((k[i]%2 == 0) and (n[i] > n_max)
                    and (graph_types[t][i] != GraphType::Chain))
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
}


/* AOF grid search: loop over N, return candidate that minimises
objective function. No need to generate complete graph
  * \param num_factors          Number of Kronecker factors
  * \param grid_center          Centre of N-k-grid
  * \param target_parameters    The target parameters
  * \param obj_func             The objective functions for each parameter
  * \param weights              The weights used in the aggregate objective
                                function
  * \param max_error            The maximum error used to determine grid
                                search bounds
  * \param log                  The model logger
  * \return result              A vector of tuples of vectors representing
                                the factor properties (N, k, type)
  * \throws invalid_argument    If target parameters are given and number of
                                target parameters does not match that of
                                objective functions or weights
*/
template<typename Logger>
std::pair<std::vector<std::tuple<factor, factor, std::vector<GraphType>>>, double> grid_search(
    const size_t& min_factors,
    const size_t& max_factors,
    const std::pair<size_t, size_t>& grid_center,
    const std::vector<double>& target_parameters,
    std::vector<std::function<double(double, factor, factor, std::vector<GraphType>)>>& obj_func, // use variadic function template with a fold expression
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
           throw invalid_argument("Number of target parameters, objective "
                                       "functions, and weights do not match");
    }

    // Get grids in N and k using interval bounds
    log->debug("Setting up grid now");
    const auto n_grid = get_N_grid(grid_center.first, max_error, min_factors, max_factors);
    const auto k_grid = get_k_grid(grid_center.second, max_error, min_factors, max_factors);
    log->debug("Grid set up");

    // Vector for resulting factors
    std::vector<std::tuple<factor, factor, std::vector<GraphType>>> res = {};

    // Initialise error
    double error = -1;

    // Grid search
    for (const auto& n: n_grid) {
         for (const auto& k: k_grid) {

            if (not check_validity(n, k)) {
                continue;
            }
            std::vector<std::vector<GraphType>> graph_types;

            // List possible type combinations
            find_possible_types(n, k, target_parameters, graph_types, max_error);

            // calculate error function for all types considered
            for (const auto& t : graph_types) {
                double current_err = 0;
                for (int i = 0; i < target_parameters.size(); ++i) {
                    current_err += weights[i]*obj_func[i](target_parameters[i], n, k, t);
                }

                // if error constant, add to set of current Pareto points
                if (current_err == error) {

                    res.emplace_back(std::make_tuple(n, k, t));

                    error = current_err;
                }
                // if error function reduced, use curent n, m, type
                else if((current_err < error) or (error == -1)) {
                    res = {{n, k, t}};
                    error = current_err;
                }
            }
        }
    }

    return {res, error};
}

/*
grid search with particle swarm
1. for each f_i, calculate l best solutions
2. ...
*/
}
#endif
