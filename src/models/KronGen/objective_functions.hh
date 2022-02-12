#ifndef UTOPIA_MODELS_KRONGEN_OBJECTIVE_FUNCS
#define UTOPIA_MODELS_KRONGEN_OBJECTIVE_FUNCS
#define UNUSED(expr) do { (void)(expr); } while (0)

#include "graph_properties.hh"
#include "type_definitions.hh"
#include "utils.hh"

namespace Utopia::Models::KronGen::ObjectiveFuncs {

using namespace Utopia::Models::KronGen;
using namespace Utopia::Models::KronGen::TypeDefinitions;


// Error norm
// TO DO: initialise this from the front end
double err_func(double x, double y) {
    return std::abs(1-x/y); // L1 norm
    //return pow((1-x/y), 2); //L2 norm
}

/// Clustering objective function for various graphs types
/**
  * \param c_t     The target clustering coefficient
  * \param graphs  A vector of graph descriptors
  * \param rng     The random number generator
  *
  * \return err    The error in the specified norm
  *
  */
double clustering_obj_func (const double& c_t,
                            const std::vector<GraphDesc>& graphs,
                            std::mt19937& rng)
{
    // Calculate graph degree variances and clustering coefficients
    double c_current=0., c_previous, k_current, k_previous, v_current, v_previous;
    size_t i = 0;
    for (const auto& graph : graphs){
        k_current = graph.mean_degree;
        c_current = GraphProperties::clustering_estimation<NWType>(
            graph.num_vertices, graph.mean_degree, graph.type, c_t, rng
        );
        v_current = GraphProperties::degree_variance(
            graph.num_vertices, graph.mean_degree, graph.type
        );
        if (i != 0) {
            c_current = Utils::Kronecker_clustering(c_current, c_previous,
                                                    k_current, k_previous,
                                                    v_current, v_previous);
            v_current = Utils::Kronecker_degree_variance(k_current, k_previous,
                                                         v_current, v_previous);
            k_current = Utils::Kronecker_mean_degree(k_current, k_previous);
        }
        else {
          i = 1;
        }
        c_previous = c_current;
        k_previous = k_current;
        v_previous = v_current;

    }

    return err_func(c_current, c_t);
}

/// Diameter objective function for various graphs types
/**
  * \param d_t     The target diameter
  * \param graphs  A vector of graph descriptors
  *
  * \return err    The error in the specified norm
  *
  */
double diameter_obj_func (const double& d_t,
                          const std::vector<GraphDesc>& graphs,
                          std::mt19937& rng)
{
    UNUSED(rng);
    double d = 0;
    for (const auto& graph : graphs){
        const size_t d_est = GraphProperties::diameter_estimation(
            graph.num_vertices, graph.mean_degree, graph.type);
        if (d_est > d) {
            d = d_est;
        }
    }

    return err_func(d, d_t);

}

/// Number of vertices objective function, independent of graph type
/**
  * \param N_t     The target vertex count
  * \param graphs  A vector of graph descriptors
  *
  * \return err    The error in the specified norm
  */
double N_obj_func (const double& N_t,
                   const std::vector<GraphDesc>& graphs,
                   std::mt19937& rng)
{
    UNUSED(rng);
    double N_res = 1;
    for (const auto& graph : graphs) {
        N_res *= graph.num_vertices;
    }

    return err_func(N_res, N_t);
}

/// Mean degree objective function for various graph types
/*
 * \param k_t     The target mean degree
 * \param graphs  A vector of graph descriptors
 *
 * \return err    The error in the specified norm
 */
double k_obj_func (const double& k_t,
                   const std::vector<GraphDesc>& graphs,
                   std::mt19937& rng)
{
    UNUSED(rng);
    double k_res = 1.0;
    for (const auto& graph : graphs) {
        k_res *= (GraphProperties::mean_degree(
          graph.num_vertices, graph.mean_degree, graph.type)+1);
    }

    return err_func(k_res, k_t+1);
}


/// Collects objective functions and weights from a config
template<typename obj_func, typename Config, typename Logger>
std::map<std::string, std::pair<double, obj_func>> get_obj_funcs(const Config& cfg,
                                                                 map_type& targets,
                                                                 const Logger& log)
{
    std::map<std::string, std::pair<double, obj_func>> objective_funcs;

    log->info("Assembling Kronecker graph with the following properties: ");
    for (const auto& p : targets) {

        const std::string param = p.first;

        // If parameter is not a target value, continue
        if (targets[param]["target"].first == false){
            continue;
        }

        const Config parameter_cfg = get_as<Config>(p.first, cfg["targets"]);

        // Collect objective functions for the various target parameters.
        // The degree sequence is not determined via an objective function
        if (param == "degree_sequence") {
            log->info("{}: {}", param, std::any_cast<std::string>(targets[param]["target"].second));
            continue;
        }
        else if (param == "clustering") {
            objective_funcs[param]
            = std::make_pair<double, obj_func>(
                get_as<double>("weight", parameter_cfg),
                clustering_obj_func);
        }
        else if (param == "diameter") {
            objective_funcs[param]
            = std::make_pair<double, obj_func>(
                get_as<double>("weight", parameter_cfg),
                diameter_obj_func);
        }
        else if (param == "mean_degree") {
            objective_funcs[param]
            = std::make_pair<double, obj_func>(
                get_as<double>("weight", parameter_cfg),
                k_obj_func);
        }
        else if (param == "num_vertices") {
            objective_funcs[param]
            = std::make_pair<double, obj_func>(
                get_as<double>("weight", parameter_cfg),
                N_obj_func);
        }
        log->info("{}: {}", param, std::any_cast<double>(targets[param]["target"].second));
    }

    return objective_funcs;
}


} // namespace KronGen::ObjectiveFuncs

#endif // UTOPIA_MODELS_KRONGEN_OBJECTIVE_FUNCS
