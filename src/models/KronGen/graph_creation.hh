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

/// Returns the Kronecker product of a list of graphs and calculates its topological
/// properties. This method does not require generating the full graph
/**
  * \tparam Graph         The graph type
  * \param cfg            A list of graphs to Kronecker together
  * \param single_factor  Whether the graph is the tensor power of a single
                          seed, and if so to what power it should be raised
  * \param build_graph    Whether the full graph needs to be built
  * \tparam RNG           The rng tpye
  * \tparam Logger        The logger type
  * \param rng            The model random number generator
  * \param log            The model logger
  * \param analysis_cfg   The analysis config, containing the list of parameters
                          to calculate
  *
  * \return Graph         The Kronecker product graphs of the factors
  */
template<typename Graph, typename RNGType, typename Logger>
Graph create_Kronecker_graph(const Config& cfg,
                             const std::pair<bool, size_t> single_factor,
                             const bool build_graph,
                             RNGType& rng,
                             const Logger& log,
                             const Config& analysis_cfg = YAML::Node())
{

    // Create initial graph and collect properties to analyse from the config
    Graph K = Utils::build_initiator_graph<Graph>(build_graph);
    auto analysis_params = Utils::get_analysis_targets(analysis_cfg);

    const auto& factor_list = get_as<Config>("Kronecker", cfg);
    const size_t num_factors = factor_list.size();

    // Generate graph as a tensor power of a single graph
    if (single_factor.first) {
        const auto& model_map = get_as<Config>("Graph1", factor_list);
        for (int i = 0; i < single_factor.second; ++i) {

            log->debug("Generating {} graph (factor {}/{}) ...",
                       get_as<std::string>("model", model_map), i+1, num_factors);
            Graph G = AuxGraphs::create_graph<Graph>(model_map, rng);

            log->debug("Calculating factor properties ... ");
            Utils::calculate_properties(G, analysis_params, log);

            // Build full graph if specified
            if (build_graph) {
                log->debug("Assembling the graph ...");
                Utils::add_self_edges(G);
                K = Utils::Kronecker_product(K, G);
            }
        }
    }

    // Generate graph from a list of different factors
    else {

        log->info("Creating Kronecker product of {} factors...", num_factors);

        for (const auto& model_map : factor_list){
            const auto factor_name = model_map.first.as<std::string>();
            const auto& factor_cfg = model_map.second;

            // Generate graph and calculate properties
            log->debug("Generating {}/{} (type: {}) ...",
                       factor_name, num_factors, get_as<std::string>("model", factor_cfg));
            Graph G = AuxGraphs::create_graph<Graph>(factor_cfg, rng);

            log->debug("Calculating factor properties ... ");
            Utils::calculate_properties(G, analysis_params, log);

            // Build full graph if specified
            if (build_graph) {
                log->debug("Assembling graph ...");
                Utils::add_self_edges(G);
                K = Utils::Kronecker_product(K, G);
            }
        }
    }

    // Write data, remove self edges, and return the graph
    Utils::write_properties(K, analysis_params);
    Utils::remove_self_edges(K);

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
  * \param cfg            A list of topological properties
  * \param build_graph    Whether to build the full graph
  * \param rng            The model rng
  * \param log            The model logger
  * \param analysis_cfg   The analysis config, containing the list of parameters to
                          calculate. Target parameters are always calculated.
  *
  * \return Graph         The resulting Kronecker graph
  */
template<typename Graph, typename RNGType, typename Logger>
Graph create_KronGen_graph(const Config& cfg,
                           const bool build_graph,
                           RNGType& rng,
                           const Logger& log,
                           const Config& analysis_cfg = YAML::Node())
{
    std::uniform_real_distribution<double> prob_distr;

    // An objective function takes a target value and a list of graph types,
    // and returns an error value
    using objective_function
    = typename std::function<double(double, factor, factor, std::vector<GraphTypes::GraphType>, RNGType&)>;

    // Get target parameters
    auto targets = Utils::get_analysis_targets(analysis_cfg, get_as<Config>("targets", cfg["KronGen"]));
    if (targets.find("num_vertices") == targets.end()) {
      throw std::invalid_argument("You must provide the target number of vertices!");
    }
    if (targets.find("mean_degree") == targets.end()) {
        targets.insert({"mean_degree",
            Utils::entry_type{{"target", std::make_pair<bool, double>(false,
                       static_cast<double>(
                       std::round(prob_distr(rng)
                       * (std::any_cast<double>(targets["num_vertices"]["target"].second)-1))))}}
        });
    }
    // Create map of objective functions and weights, and output info message.
    // All weights are 1.0 for now.
    std::map<std::string, std::pair<double, objective_function>> objective_funcs;
    log->info("Assembling Kronecker graph with the following properties: ");
    for (const auto& p : targets) {

        const std::string param = p.first;
        // If parameter is not a target value, continue
        if (targets[param]["target"].first == false){
            continue;
        }

        // Degree sequence is not determined via an objective function
        if (param == "degree_sequence") {
            log->info("{}: {}", param, std::any_cast<std::string>(targets[param]["target"].second));
            continue;
        }
        // Collect objective functions for the various target parameters
        if (param == "clustering") {
            objective_funcs[param]
            = std::make_pair<double, objective_function>(1.0, ObjectiveFuncs::clustering_obj_func);
        }
        else if (param == "diameter") {
            objective_funcs[param]
            = std::make_pair<double, objective_function>(1.0, ObjectiveFuncs::diameter_obj_func);
        }
        else if (param == "mean_degree") {
            objective_funcs[param]
            = std::make_pair<double, objective_function>(1.0, ObjectiveFuncs::k_obj_func);
        }
        else if (param == "num_vertices") {
            objective_funcs[param]
            = std::make_pair<double, objective_function>(1.0, ObjectiveFuncs::N_obj_func);
        }
        log->info("{}: {}", param, std::any_cast<double>(targets[param]["target"].second));
    }

    // Get grid search parameters: get minimum and maximum factor size
    const size_t d_min = get_as<size_t>("min_dimension", cfg["KronGen"]);
    const size_t d_max = get_as<size_t>("max_dimension", cfg["KronGen"]);

    // Get grid center
    const auto N_0 = std::any_cast<double>(targets["num_vertices"]["target"].second);
    const auto k_0 = std::any_cast<double>(targets["mean_degree"]["target"].second);
    const auto grid_center = std::make_pair<size_t, size_t>(N_0, k_0);

    // Get grid search scope
    const double tolerance = get_as<double>("tolerance", cfg["KronGen"]);

    // Perform grid search
    const auto Paretos = GridSearch::grid_search(d_min, d_max, grid_center,
                                             tolerance, targets, objective_funcs,
                                             log, rng);

    // Collect the number of Pareto points found
    const size_t n_Paretos = Paretos.size();
    targets["num_Paretos"]["calculate"].second = n_Paretos;
    if (n_Paretos == 0) {
        log->warn("No Pareto points found. It may be that your configuration "
                  "settings are invalid. Consider changing the target parameters"
                  " or grid size."
         );
         return Utils::build_initiator_graph<Graph>(false);
    }
    log->info("Grid search complete. Found {} Pareto point(s).", n_Paretos);

    // If multiple Pareto points found, randomly pick one
    const size_t rdm_pt = std::round(prob_distr(rng) * (n_Paretos-1));
    const auto N_factors = std::get<0>(Paretos[rdm_pt]);
    const auto k_factors = std::get<1>(Paretos[rdm_pt]);
    const auto types = std::get<2>(Paretos[rdm_pt]);

    // Keep track of number of non-trivial factors
    auto n_factors = N_factors.size();
    targets["num_factors"]["calculate"].second = n_factors;
    log->info("Number of Kronecker factors: {}.", n_factors);

    // ... Create initiator graph
    Graph K = Utils::build_initiator_graph<Graph>(build_graph);

    // Assemble Kronecker graph
    for (size_t i = 0; i < N_factors.size(); ++i) {

        // Get graph generation parameters
        const auto n_curr = N_factors[i];
        const auto k_curr = k_factors[i];

        // Ignore trivial factors
        if (k_curr == 0 or n_curr == 1) {
            targets["num_factors"]["calculate"].second = --n_factors;
            continue;
        }

        // Keep track of largest factor size
        if (n_curr > std::any_cast<size_t>(targets["largest_comp"]["calculate"].second)) {
            targets["largest_comp"]["calculate"].second = n_curr;
        }

        // Create graphs
        double c_t = std::any_cast<double>(targets["clustering"]["target"].second);
        Graph H = AuxGraphs::create_graph<Graph>(n_curr, k_curr, types[i], rng, c_t);

        log->debug("Current factor: {} graph with {} vertices, mean degree k = {}",
                   GraphProperties::Graph_Type[types[i]], n_curr, k_curr);

        // Calculate properties of Kronecker product
        Utils::calculate_properties(H, targets, log);

        // Build full graph if specified
        if (build_graph) {
            log->debug("Assembling graph ...");
            Utils::add_self_edges(H);
            K = Utils::Kronecker_product(K, H);
        }
    }

    // Write data, remove self edges, and return the graph
    log->debug("Writing properties ...");
    Utils::write_properties(K, targets);
    Utils::remove_self_edges(K);

    log->info("Graph creation complete.");

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
    Config nw_cfg = YAML::Node();
    try {
        nw_cfg = get_as<Config>("graph_analysis", cfg["NetworkAnalyser"]);
        includes_analysis_cfg = true;
    }
    catch (YAML::InvalidNode&){}
    catch (Utopia::KeyError&){}
    if (includes_analysis_cfg) {
        const bool analysis_enabled = get_as<bool>("enabled", nw_cfg);
        if (not analysis_enabled) {
            nw_cfg = YAML::Node();
        }
    }

    // Get graph creation config
    const Config graph_cfg = includes_analysis_cfg
        ? get_as<Config>("create_graph", cfg)
        : cfg;

    // Get the graph generating model
    const std::string model = get_as<std::string>("model", graph_cfg);

    if (model == "Kronecker") {
        const bool build_graph = get_as<bool>("build_graph", graph_cfg);
        const size_t num_factors = get_as<Config>("Kronecker", graph_cfg).size();
        std::pair<bool, size_t> single_factor{false, 0};
        try {
            single_factor = {true, get_as<size_t>("power", graph_cfg["Kronecker"])};
        }
        catch (Utopia::KeyError&){}

        if (single_factor.first and (num_factors > 2)){
            throw std::invalid_argument("You have specified both a tensor power"
                      " graph and multiple factor graphs! Either remove one of "
                      " the graphs or remove the 'power' argument!");
        }

        return create_Kronecker_graph<Graph>(graph_cfg, single_factor, build_graph,
                                             rng, log, nw_cfg);
    }

    else if (model == "KronGen") {
        const bool build_graph = get_as<bool>("build_graph", graph_cfg);
        return create_KronGen_graph<Graph>(graph_cfg, build_graph, rng, log, nw_cfg);
    }

    // If the graph is not a Kronecker product graph, directly generate the graph
    // using a standard graph generation algorithm. The AuxGraphs::create_graph
    // function contains all the Utopia graph generation algorithms and more
    else {
        return AuxGraphs::create_graph<Graph>(graph_cfg, rng);
    }
}

} // namespace KronGen::GraphCreation

#endif // UTOPIA_MODELS_KRONGEN_GRAPHCREATION
