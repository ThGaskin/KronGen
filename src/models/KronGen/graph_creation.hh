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
#include "clustering.hh"
#include "diameter.hh"
#include "../NetworkAnalyser/graph_metrics.hh"
#include "utils.hh"

namespace Utopia::Models::KronGen::GraphCreation {

using namespace boost;
using namespace Utopia::Models::KronGen;
using namespace Utopia::Models::NetworkAnalyser;

/// Returns the Kronecker product of a list of graphs
/**
  * \param cfg            A list of graphs to Kronecker together
  * \param rng            The model rng
  * \param log            The model logger
  * \param analysis_cfg   The analysis config, containing the list of parameters to
                          calculate
*/
template<typename Graph, typename RNGType, typename Logger>
Graph create_Kronecker_graph(const Config& cfg,
                             RNGType& rng,
                             const Logger& log,
                             const Config& analysis_cfg = YAML::Node(YAML::NodeType::Map))
{
    log->info("Creating Kronecker graph");
    std::uniform_real_distribution<double> distr(0, 1);

    // ... Create graph with one vertex and a self-edge ........................
    Graph K{1};
    Utils::add_self_edges(K);

    // ... Data containers for graph analysis properties .......................
    double c_global = -1;
    double diam = -1;
    double mean_deg;
    double variance;

    // ... Get analysis targets ................................................
    bool calculate_c;
    try {
        calculate_c = (get_as<bool>("clustering_global", analysis_cfg)
                    && get_as<bool>("enabled", analysis_cfg));
    }
    catch (YAML::InvalidNode&){}
    catch (Utopia::KeyError&){}

    bool calculate_diam;
    try {
        calculate_diam = (get_as<bool>("diameter", analysis_cfg)
                    && get_as<bool>("enabled", analysis_cfg));
    }
    catch (YAML::InvalidNode&){}
    catch (Utopia::KeyError&){}

    bool first_run = true;

    // Generate Graph
    for (const auto& model_map : cfg["Kronecker"]){
        const auto& model_cfg = model_map.second;
        Graph G = AuxGraphs::create_graph<Graph>(model_cfg, rng);

        // ... Calulate properties: stored in first vertex .....................
        Utils::calculate_properties (G,
                                     first_run,
                                     calculate_c,
                                     calculate_diam,
                                     c_global,
                                     diam,
                                     mean_deg,
                                     variance);
        // .....................................................................

        Utils::add_self_edges(G);

        K = Utils::Kronecker_product(K, G, rng, distr);
    }

    // Remove self-loops
    Utils::remove_self_edges(K);

    // Write data
    if (calculate_c) {
        K[0].state.clustering_global = c_global;
    }
    if(calculate_diam) {
        K[0].state.diameter = diam;
    }

    // Return the graph
    return K;

}

/// Creates a graph from a list of topological properties
/**
  * \param cfg            A list of topological properties and tolerances
  * \param rng            The model rng
  * \param log            The model logger
  * \param analysis_cfg   The analysis config, containing the list of parameters to
                          calculate
*/
template<typename Graph, typename RNGType, typename Logger>
Graph create_KronGen_graph(const Config& cfg,
                           RNGType& rng,
                           const Logger& log,
                           const Config& analysis_cfg = YAML::Node(YAML::NodeType::Map))
{
    std::uniform_real_distribution<double> distr(0, 1);

    // ... Get topological properties ..........................................
    const double N = get_as<double>("num_vertices", cfg);
    const double m = get_as<double>("mean_degree", cfg);
    double c = -1;
    std::string degree_distr = "";
    double diameter = -1;

    try {
        c = get_as<double>("clustering_coeff", cfg["KronGen"]);
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

    // ... Get analysis targets ................................................
    bool calculate_c = false;
    if (c != -1) {
        calculate_c = true;
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
    }
    else {
        try {
            calculate_diam = (get_as<bool>("diameter", analysis_cfg)
                        && get_as<bool>("enabled", analysis_cfg));
        }
        catch (YAML::InvalidNode&){}
        catch (Utopia::KeyError&){}
    }

    // ... Graph creation ......................................................
    // Output info message with given properties
    log->info("Assembling KronGen {} graph with {} vertices, mean degree m = {}"
               "{}{} ...", degree_distr, N, m,
               (c != -1 ? ", clustering coefficient = "+to_string(c) : ""),
               (diameter != -1 ? ", diameter = "+to_string(diameter) : ""));

    // ... Create graphs when no properties are passed .........................
    if ((diameter == -1) and (c == -1)) {

        if (degree_distr == "scale-free") {
            return Utopia::Graph::create_BarabasiAlbert_graph<Graph>(N, m,
                                                                     false, rng);
        }
        else {
            return Utopia::Graph::create_ErdosRenyi_graph<Graph>(N, m, false,
                                                                 false, rng);
        }
    }

    // ... Create graph with given diameter ....................................
    Graph K{1};
    Utils::add_self_edges(K);

    if (diameter > 0) {

        Diameter::create_diameter_graph(K, N, m, c, diameter, degree_distr,
                                        calculate_c, calculate_diam,
                                        rng, distr, log);
    }

    // ... Create graph with given clustering coefficient ......................
    if (c != -1) {
        const double tolerance = get_as<double>("tolerance", cfg["KronGen"]);
        Clustering::create_clustering_graph(K, N, m, c, diameter, degree_distr,
                                            calculate_c, calculate_diam, tolerance,
                                            rng, distr, log);
    }

    // .........................................................................

    Utils::remove_self_edges(K);

    // Return the graph
    return K;

}

// .............................................................................

/// Custom create_graph function
/**
  * \param cfg            Graph config
  * \param rng            The model rng
  * \param log            The model logger
  * \param analysis_cfg   The analysis config, containing the list of parameters to
                          calculate
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
        return create_Kronecker_graph<Graph>(graph_cfg, rng, log, nw_cfg);
    }

    else if (model == "KronGen") {
        return create_KronGen_graph<Graph>(graph_cfg, rng, log, nw_cfg);
    }

    else {
        return AuxGraphs::create_graph<Graph>(graph_cfg, rng);
    }

}

}
#endif
