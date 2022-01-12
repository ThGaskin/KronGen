#ifndef UTOPIA_MODELS_KRONGEN_GRAPHCREATION
#define UTOPIA_MODELS_KRONGEN_GRAPHCREATION

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
#include "graph_properties.hh"
#include "graph_types.hh"
#include "grid_search.hh"
#include "objective_functions.hh"
#include "../NetworkAnalyser/graph_metrics.hh"
#include "utils.hh"

namespace Utopia::Models::KronGen::GraphCreation {

using namespace Utopia::Models::KronGen;

using factor = typename std::vector<std::size_t>;
using factors = typename std::vector<std::vector<std::size_t>>;

/// Returns the Kronecker product of a list of graphs
/**
  * \tparam Graph         The graph type
  * \tparam RNG           The rng tpye
  * \tparam Logger        The logger type
  *
  * \param cfg            A list of graphs to Kronecker together
  * \param rng            The model random number generator
  * \param log            The model logger
  * \param analysis_cfg   The analysis config, containing the list of parameters
                          to calculate
  *
  * \return Graph         The Kronecker product graphs of the factors
  */
template<typename Graph, typename RNGType, typename Logger>
Graph create_Kronecker_graph(const Config& cfg,
                             const int tensor_power,
                             RNGType& rng,
                             const Logger& log,
                             const Config& analysis_cfg = YAML::Node(YAML::NodeType::Map))
{
    log->info("Creating Kronecker graph ...");

    // ... Create initator graph ...............................................
    const bool build_graph = get_as<bool>("build_graph", cfg);
    Graph K = Utils::build_initiator_graph<Graph>(build_graph);

    // ... Data containers for graph analysis properties .......................
    double c = -1;
    double diam = -1;
    double mean_deg;
    size_t num_vertices;
    double variance;

    // ... Get analysis targets ................................................
    bool calculate_c = false;
    try {
        calculate_c = (get_as<bool>("clustering_global", analysis_cfg)
                    && get_as<bool>("enabled", analysis_cfg));
    }
    catch (YAML::InvalidNode&){}
    catch (Utopia::KeyError&){}

    bool calculate_diam = false;
    try {
        calculate_diam = (get_as<bool>("diameter", analysis_cfg)
                    && get_as<bool>("enabled", analysis_cfg));
    }
    catch (YAML::InvalidNode&){}
    catch (Utopia::KeyError&){}

    bool first_run = true;

    const size_t num_factors = get_as<Config>("Kronecker", cfg).size();
    if (tensor_power > 1 and (num_factors > 2)){
        throw std::invalid_argument("You have specified both a tensor power graph"
                  " and multiple graph types! Either remove one of the graphs"
                  " or remove the 'power' argument!");
    }
    // Generate Graph
    if (tensor_power < 0) {
        int i = 0;
        for (const auto& model_map : cfg["Kronecker"]){
            const auto& model_cfg = model_map.second;
            log->debug("Generating {} graph (factor {}/{}) ...",
                       get_as<std::string>("model", model_cfg), i+1, num_factors);
            Graph G = AuxGraphs::create_graph<Graph>(model_cfg, rng);

            // ... Calulate properties: stored in first vertex .................
            Utils::calculate_properties (G,
                                         first_run,
                                         calculate_c,
                                         calculate_diam,
                                         c,
                                         diam,
                                         mean_deg,
                                         num_vertices,
                                         variance);
            // .................................................................
            if (build_graph) {
                Utils::add_self_edges(G);
                K = Utils::Kronecker_product(K, G);
            }
            ++i;
        }
    }
    else {
        const auto& model_cfg = get_as<Config>("Kronecker", cfg);
        const auto& model_map = get_as<Config>("Graph1", model_cfg);
        for (int i = 0; i<tensor_power; ++i) {
            log->debug("Generating {} graph (factor {}/{}) ...",
                       get_as<std::string>("model", model_map), i+1, tensor_power);
            Graph G = AuxGraphs::create_graph<Graph>(model_map, rng);

            // ... Calulate properties: stored in first vertex .................
            Utils::calculate_properties (G,
                                         first_run,
                                         calculate_c,
                                         calculate_diam,
                                         c,
                                         diam,
                                         mean_deg,
                                         num_vertices,
                                         variance);
            // .................................................................
            if (build_graph) {
                Utils::add_self_edges(G);
                K = Utils::Kronecker_product(K, G);
            }
        }
    }

    // Remove self-loops
    Utils::remove_self_edges(K);

    // Write data
    if (calculate_c) {
        K[0].state.clustering_global = c;
    }
    if (calculate_diam) {
        K[0].state.diameter = diam;
    }
    K[0].state.mean_deg = mean_deg;
    K[0].state.num_vertices = num_vertices;

    // Return the graph
    return K;

}

/// Creates a Kronecker graph from a list of topological properties by conducting
/// a grid search over the specified ensemble, optimising an aggregate objective
/// function, and returning the Pareto front.
/**
  * \tparam Graph         The graph type
  * \tparam RNG           The rng tpye
  * \tparam Logger        The logger type
  *
  * \param cfg            A list of topological properties and tolerances
  * \param rng            The model rng
  * \param log            The model logger
  * \param analysis_cfg   The analysis config, containing the list of parameters to
                          calculate
  *
  * \return Graph         The resulting Kronecker graph
  */
template<typename Graph, typename RNGType, typename Logger>
Graph create_KronGen_graph(const Config& cfg,
                           RNGType& rng,
                           const Logger& log,
                           const Config& analysis_cfg = YAML::Node(YAML::NodeType::Map))
{

    // ... Get topological properties ..........................................
    const double N = get_as<double>("num_vertices", cfg);
    const double k = get_as<double>("mean_degree", cfg);
    double c = -1;
    std::string degree_distr = "";
    double diameter = -1;

    try {
        c = get_as<double>("clustering", cfg["KronGen"]);
        if (c < 0 or c > 1) {
            throw std::invalid_argument("Clustering coefficient must be in [0, 1]!");
        }
    }
    catch(Utopia::KeyError&){}

    try {
        degree_distr = get_as<std::string>("degree_distribution", cfg["KronGen"]);
    }
    catch(Utopia::KeyError&){}

    try {
        diameter = get_as<double>("diameter", cfg["KronGen"]);
        if (diameter <= 0) {
            throw std::invalid_argument("Diameter must be greater than 0!");
        }
    }
    catch(Utopia::KeyError&){}

    // ... Collect target parameters, objective functions, and weights .........
    std::vector<double> target_params = {N, k};
    std::vector<std::function<double(double, factor, factor,
        std::vector<GraphTypes::GraphType>)>> obj_func =
                                             {ObjectiveFuncs::N_obj_func,
                                              ObjectiveFuncs::k_obj_func};
    std::vector<double> weights = {1, 1};

    bool calculate_c = false;
    if (c != -1) {
        calculate_c = true;
        obj_func.emplace_back(ObjectiveFuncs::clustering_obj_func);
        target_params.emplace_back(c);
        weights.emplace_back(1.);
    }
    else {
        try {
            calculate_c = (get_as<bool>("clustering_global", analysis_cfg)
                        && get_as<bool>("enabled", analysis_cfg));
        }
        catch (YAML::InvalidNode&){}
        catch (Utopia::KeyError&){}
    }

    bool calculate_diam = false;
    if (diameter != -1) {
        calculate_diam = true;
        obj_func.emplace_back(ObjectiveFuncs::diameter_obj_func);
        target_params.emplace_back(diameter);
        weights.emplace_back(1.);
    }
    else {
        try {
            calculate_diam = (get_as<bool>("diameter", analysis_cfg)
                        && get_as<bool>("enabled", analysis_cfg));
        }
        catch (YAML::InvalidNode&){}
        catch (Utopia::KeyError&){}
    }

    // Output info message with given properties
    log->info("Assembling KronGen {} graph with {} vertices, mean degree k = {}"
               "{}{} ...", degree_distr, N, k,
               (c != -1 ? ", clustering coefficient = "+std::to_string(c) : ""),
               (diameter != -1 ? ", diameter = "+std::to_string(diameter) : ""));

    // ... Create graphs when no properties are passed .........................
    if ((diameter == -1) and (c == -1)) {

        if (degree_distr == "scale-free") {
            return Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N, k,
                                                                     false, rng);
        }
        else {
            return Utopia::Graph::create_ErdosRenyi_graph<Graph>(N, k, false,
                                                                 false, rng);
        }
    }

    // Ignore graphs with specified degree distribution for now
    else if (degree_distr != "") {

      return Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N, k, false, rng);

    }

    // ... Create graph with given properties via grid search ..................
    // Get grid search parameters
    const size_t d_min = get_as<size_t>("min_dimension", cfg["KronGen"]);
    const size_t d_max = get_as<size_t>("max_dimension", cfg["KronGen"]);
    const double tolerance = get_as<double>("tolerance", cfg["KronGen"]);
    const auto grid_center = std::make_pair(N, k);

    // Perform grid search
    const auto res = GridSearch::grid_search(d_min,
                                             d_max,
                                             grid_center,
                                             target_params,
                                             obj_func,
                                             weights,
                                             tolerance,
                                             log);

    // Collect the resulting Pareto points (consisting of N, k, and t vectors)
    const auto Paretos = res.first;

    // Collect the resulting error value (equal for all Pareto points)
    const double error = res.second;

    const size_t n_Paretos = Paretos.size();
    log->info("Grid search complete. Found {} Pareto point(s).", n_Paretos);

    // If multiple Pareto points found, randomly pick one
    std::uniform_real_distribution<double> prob_distr;
    const size_t rdm_pt = std::round(prob_distr(rng) * (n_Paretos-1));
    const auto N_factors = std::get<0>(Paretos[rdm_pt]);
    const auto k_factors = std::get<1>(Paretos[rdm_pt]);
    const auto types = std::get<2>(Paretos[rdm_pt]);

    // Keep track of number of non-trivial factors
    auto n_factors = N_factors.size();

    log->info("Number of Kronecker factors: {}.", n_factors);

    // ... Create initiator graph
    const bool build_graph = get_as<bool>("build_graph", cfg["KronGen"]);
    Graph K = Utils::build_initiator_graph<Graph>(build_graph);

    // Keep track of resulting properties as Kronecker graph is assembled
    double c_res = -1;
    double diam_res = -1;
    double k_res = 0;
    double N_res = 1;
    double var_res = -1;
    size_t largest_factor = 0;

    for (size_t i = 0; i < N_factors.size(); ++i) {

        // Get graph generation parameters
        const auto n_curr = N_factors[i];

        // Calculate properties
        double k_curr = GraphProperties::mean_degree(n_curr, k_factors[i], types[i]);
        double c_curr = GraphProperties::clustering_estimation(n_curr, k_curr, types[i]);
        double diam_curr = GraphProperties::diameter_estimation(n_curr, k_curr, types[i]);
        double var_curr = GraphProperties::degree_variance(n_curr, k_curr, types[i]);

        // Ignore trivial factors
        if (k_curr == 0 or n_curr == 1) {
          --n_factors;
          continue;
        }
        if (n_curr > largest_factor) {
            largest_factor = n_curr;
        }

        // Create graphs and correct estimated properties if necessary
        Graph H;

        if (types[i] == GraphTypes::Chain) {
            H = AuxGraphs::create_chain_graph<Graph>(n_curr);
        }
        else if (types[i] == GraphTypes::Complete){
            H = Utopia::Graph::create_complete_graph<Graph>(n_curr);
        }
        else if (types[i] == GraphTypes::ErdosRenyi) {
            H = Utopia::Graph::create_ErdosRenyi_graph<Graph>(n_curr, k_curr,
                                                         false, false, rng);
            // hm ...
            if (k_curr < 10) {
                if (calculate_c) {
                    c_curr = Utopia::Models::NetworkAnalyser::global_clustering_coeff(H);
                }
                var_curr = Utopia::Models::NetworkAnalyser::degree_statistics(H).second;
            }

            if (calculate_diam) {
                diam_curr = Utopia::Models::NetworkAnalyser::diameter(H);
            }

        }
        else if (types[i] == GraphTypes::Regular) {
            H = Utopia::Graph::create_regular_graph<Graph>(n_curr, k_curr, false);
            if (calculate_diam) {
                diam_curr = Utopia::Models::NetworkAnalyser::diameter(H);
            }
        }

        log->debug("Current factor: {} graph with {} vertices, mean degree k = {}"
                   "{}{} ...", GraphProperties::Graph_Type[types[i]], n_curr, k_curr,
                   (c != -1 ? ", clustering coefficient = "+std::to_string(c_curr) : ""),
                   (diameter != -1 ? ", diameter = "+std::to_string(diam_curr) : ""));

        // ... Calculate properties of Kronecker product .......................
        if (calculate_c) {
            if (c_res == -1) {
                c_res = c_curr;
                var_res = var_curr;
            }
            else {
                c_res = Utils::Kronecker_clustering(c_res, c_curr,
                                                    k_res, k_curr,
                                                    var_res, var_curr);
                var_res = Utils::Kronecker_degree_variance(k_res, k_curr,
                                                           var_res, var_curr);

            }
        }
        if (calculate_diam) {
            diam_res = std::max(diam_res, diam_curr);
        }
        k_res = Utils::Kronecker_mean_degree(k_res, k_curr);
        N_res = Utils::Kronecker_num_vertices(N_res, n_curr);

        // If full graph being generated: add self-edges and create Kronecker product
        if (build_graph) {
            Utils::add_self_edges(H);
            K = Utils::Kronecker_product(K, H);
        }
    }

    // ... Write properties ....................................................
    K[0].state.clustering_global = c_res;
    K[0].state.diameter = diam_res;
    K[0].state.error = error;
    K[0].state.largest_comp = largest_factor;
    K[0].state.mean_deg = k_res;
    K[0].state.num_factors = n_factors;
    K[0].state.num_Paretos = n_Paretos;
    K[0].state.num_vertices = N_res;
    K[0].state.var = var_res;

    Utils::remove_self_edges(K);

    // Return the graph
    return K;

}

// .. Custom graph creation function ...........................................

/// Create a graph from a configuration node
/** Select a graph creation algorithm and create the graph object a
 *  configuration node.
 *
 * \tparam Graph        The graph type
 * \tparam RNG          The random number generator type
 * \tparam Logger       The logger type
 *
 * \param cfg           The configuration
 * \param rng           The random number generator
 * \param log           The model logger
 * \param includes_analysis_cfg   Whether or not the config includes configuration
                        which lists all properties to analyse
 *
 * \return Graph        The graph
 */
template<typename Graph, typename RNGType, typename Logger>
Graph create_graph(const Config& cfg,
                   RNGType& rng,
                   const Logger& log,
                   bool includes_analysis_cfg = false)
{

    // Get graph analysis config, if provided
    Config nw_cfg = YAML::Node(YAML::NodeType::Map);
    try {
        nw_cfg = get_as<Config>("graph_analysis", cfg["NetworkAnalyser"]);
        includes_analysis_cfg = true;
    }
    catch (YAML::InvalidNode&){}
    catch (Utopia::KeyError&){}

    // Get graph creation config
    const Config graph_cfg = includes_analysis_cfg
        ? get_as<Config>("create_graph", cfg)
        : cfg;

    // Get the graph generating model
    const std::string model = get_as<std::string>("model", graph_cfg);

    if (model == "Kronecker") {
        int tensor_power = -1;
        try {
            tensor_power = get_as<int>("power", graph_cfg["Kronecker"]);
        }
        catch (YAML::InvalidNode&){}
        catch (Utopia::KeyError&){}
        return create_Kronecker_graph<Graph>(graph_cfg, tensor_power, rng, log, nw_cfg);
    }

    else if (model == "KronGen") {
        return create_KronGen_graph<Graph>(graph_cfg, rng, log, nw_cfg);
    }

    // If the graph is not a Kronecker product graph, directly generate graph
    // using a graph generation algorithm
    else {
        return AuxGraphs::create_graph<Graph>(graph_cfg, rng);
    }
}

} // namespace KronGen::GraphCreation

#endif // UTOPIA_MODELS_KRONGEN_GRAPHCREATION
