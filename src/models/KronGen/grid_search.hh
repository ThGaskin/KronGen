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
#include "utils.hh"

#include "../NetworkAnalyser/graph_metrics.hh"

namespace Utopia::Models::KronGen::GridSearch {

using namespace std;
using namespace boost;
using namespace Utopia::Models::KronGen;
using namespace Utopia::Models::KronGen::Utils;
using namespace Utopia::Models::NetworkAnalyser;

using factor = typename std::vector<std::size_t>;

/* AOF grid search for ER graphs: loop over N, return candidate that minimises
objective function. No need to generate complete graph
  * \param num_factors          Number of Kronecker factors
  * \param grid_center          Centre of N-k-grid
  * \param target_parameters    The target parameters
  * \param obj_func             The objective functions for each parameter
  * \param weights              The weights used in the aggregate objective
                                function
  * \param max_err              The maximum error used to determine grid
                                search bounds
  * \return result              A vector of pairs of vectors representing
                                the factor properties
  * \throws invalid_argument    If target parameters are given and number of
                                target parameters does not match that of
                                objective functions or weights
*/
template<typename Logger>
std::vector<std::pair<factor, factor>> grid_search_ER(
    const size_t& min_factors,
    const size_t& max_factors,
    const std::pair<size_t, size_t>& grid_center,
    const std::vector<double>& target_parameters,
    std::vector<std::function<double(double, factor, factor)>>& obj_func, // use variadic function template with a fold expression
    const std::vector<double>& weights,
    const double& max_error,
    const Logger& log,
    const std::pair<int, int>& diameter_params,
    std::vector<std::pair<bool, int>>& chain_props)
{

    // If no target parameters are given, simply return starting values
    if (target_parameters.empty()) {
        return {{{grid_center.first}, {grid_center.second}}};
    }

    std::vector<std::pair<factor, factor>> res = {};

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

    // Initialise error
    double error = -1;

    // Grid search
    for (const auto& n: n_grid) {
         for (const auto& k: k_grid) {

            if (n.size() != k.size()) {
                continue;
            }
            // if target parameter is diameter and n[i] == d
            // use chain graph

            bool invalid_config = false;
            bool chain = false;
            std::pair<bool, int> chain_property = {false, -1};
            // mean degrees must be smaller than number of vertices
            for (int i = 0; i < k.size(); ++i) {

                if (k[i]>=n[i]){
                    invalid_config = true;
                    break;
                }
                if (k[i]<=1 and n[i]!=k[i]+1) {
                    invalid_config = true;
                    break;
                }

                // ... chain ...................................................
                // should allow chains of any length (and any number of chains?)
                // the chain has clustering zero - need a correction term in the error function
                // the chain does not have mean degree 2, but rather 2(N-1)/N - need a correction in the error function
                // if chain does not change error - needs to be added to pareto set
                if ((n[i] == diameter_params.second+1) and (k[i]==2) and (not chain)) {
                    chain = true;
                    chain_property = {true, i};
                }
                // .............................................................


            }
            if (invalid_config){continue;}

            // calculate error function
            double current_err = 0;
            for (int i = 0; i < target_parameters.size(); ++i) {

                // ... chain ...................................................
                if (i == diameter_params.first and chain) {
                    continue; // problem: chain changes the error of clustering and mean degree!
                }
                // .............................................................

                current_err += weights[i]*obj_func[i](target_parameters[i], n, k);
            }
            // if error constant, add to set of current Pareto points
            if (current_err == error) {

                res.emplace_back(std::make_pair(n, k));

                // .............................................................
                if (diameter_params.first != -1) {
                    chain_props.emplace_back(chain_property);
                }
                // .............................................................

                error = current_err;
            }
            // if error function reduced, use curent n, m
            else if((current_err < error) or (error == -1)) {
                res = {{n, k}};

                // .............................................................
                if (diameter_params.first != -1) {
                    chain_props = {chain_property};
                }
                // .............................................................

                error = current_err;
            }
        }
    }

    return res;
}

/*
grid search with particle swarm
1. for each f_i, calculate l best solutions
2. ...
*/
}
#endif
