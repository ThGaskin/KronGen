#ifndef UTOPIA_MODELS_KRONGEN_UTILS
#define UTOPIA_MODELS_KRONGEN_UTILS

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include "utopia/core/graph/iterator.hh"

#include "type_definitions.hh"
#include "../NetworkAnalyser/graph_metrics.hh"

namespace Utopia::Models::KronGen::Utils {

using namespace Utopia::Models::KronGen::TypeDefinitions;

// .............................................................................
// ... Kronecker product utility functions .....................................
// .............................................................................

/// Kronecker product of undirected graphs. Graphs must have a self-loop on every node
/**
  * \tparam  Graph  The graph type
  * \param G        The first factor
  * \param H        The second factor
  *
  * \return K       The Kronecker product
  */
template<typename Graph>
Graph Kronecker_product(Graph& G, Graph& H)
{
    const size_t M = boost::num_vertices(H);
    Graph K{boost::num_vertices(G) * M};

    for(const auto g : Utopia::range<IterateOver::edges>(G)) {
        for(const auto h : Utopia::range<IterateOver::edges>(H)) {
            const auto s_1 = boost::source(g, G);
            const auto s_2 = boost::source(h, H);
            const auto t_1 = boost::target(g, G);
            const auto t_2 = boost::target(h, H);

            auto i = s_1*M + s_2;
            auto j = t_1*M + t_2;
            if ((not boost::edge(i, j, K).second)
                 and (not boost::edge(j, i, K).second)){
                boost::add_edge(i, j, K);
            }

            i = s_1*M + t_2;
            j = t_1*M + s_2;
            if ((not boost::edge(i, j, K).second)
                 and (not boost::edge(j, i, K).second)){
                boost::add_edge(i, j, K);
            }
        }
    }

    return K;
}

/// Calculate the number of vertices of a Kronecker product of two graphs G, H
/**
  * \param N    The number of vertices of G
  * \param M    The number of vertices of H
  *
  * \return     The number of vertices of G (x) H
*/
double Kronecker_num_vertices (const double N, const double M)
{
    return N * M;
}

/// Calculate the number of vertices of a Kronecker product of two graphs G, H
/**
  * \param g    The mean degree of G
  * \param h    The mean degree of H
  *
  * \return     The mean_degree of G (x) H
*/
double Kronecker_mean_degree (const double g, const double h)
{
    return ((g+1)*(h+1) - 1);
}

/// Calculate the degree distribution variance of a Kronecker product of two graphs G, H
/**
  * \param g    The mean degree of G
  * \param h    The mean degree of H
  * \param v_g  The degree variance of G
  * |param v_h  The degree variance of H
  *
  * \return     The degree variance of G (x) H
*/
double Kronecker_degree_variance (const double g,
                                  const double h,
                                  const double v_g,
                                  const double v_h)
{
    return (v_g*v_h + pow(1+g, 2)*v_h + pow(1+h, 2)*v_g);
}

/// Calculate the degree sequence of a Kronecker product of two graphs G, H.
/// A degree sequence is a vector of pairs {{k, n_k}}
/**
  * \param seq_G   The degree sequence of G
  * \param seq_H   The degree sequence of H
  *
  * \return       The degree sequence of G (x) H
*/
vector_pt Kronecker_degree_sequence(const vector_pt& seq_G,
                                    const vector_pt& seq_H)
{
    const size_t k_max = (seq_G.back().first+1)*(seq_H.back().first+1);
    std::vector<size_t> indices(k_max, 0);
    for (size_t i = 0; i < seq_G.size(); ++i){
        for (size_t j = 0; j < seq_H.size(); ++j){
            const auto k = (seq_G[i].first+1)*(seq_H[j].first+1)-1;
            indices[k] += seq_G[i].second*seq_H[j].second;
        }
    }
    vector_pt seq_K;
    for (size_t i = 0; i < indices.size(); ++i){
        if (indices[i] > 0) {
            seq_K.push_back({i, indices[i]});
        }
    }

    return seq_K;
}

/// Calculate the clustering coefficient of a Kronecker product of two graphs G, H
/**
  * \param c_G        The clustering coefficient of G
  * \param c_H        The clustering coefficient of H
  * \param g          The mean degree of G
  * \param h          The mean degree of H
  * \param var_g      The degree distribution variance of G
  * \param var_h      The degree distribution variance of H
  *
  * \return           The clustering coefficient of the Kronecker product
  */
double Kronecker_clustering (const double c_G,
                             const double c_H,
                             const double g,
                             const double h,
                             const double v_g,
                             const double v_h)
{
    // Undefined clustering coefficient only occurs if a single vertex is
    // combined with a double-handle graph, both of which have undefined
    // clustering coefficients
    if (c_G == 0 and c_H == 0 and v_g == 0 and v_h == 0
        and ((g == 0 and h == 1) or (g == 1 and h == 0))){
        return 0;
    }

    const double k = Kronecker_mean_degree(g, h);
    const double v_k = Kronecker_degree_variance(g, h, v_g, v_h);

    const double t_1 = c_G*(v_g+g*(g-1)) + 3*g + 1;
    const double t_2 = c_H*(v_h+h*(h-1)) + 3*h + 1;
    const double t_3 = 3*k + 1;
    const double t_4 = (v_k) + k*(k-1);

    return ((t_1*t_2 - t_3) / t_4);
}

/// Add self_edges to a graph
template<typename Graph>
void add_self_edges (Graph& G) {
    for (const auto& v : range<IterateOver::vertices>(G)) {
        boost::add_edge(v, v, G);
    }
}

/// Remove self_edges from a graph
template<typename Graph>
void remove_self_edges (Graph& G) {
    for (const auto& v : range<IterateOver::vertices>(G)) {
        boost::remove_edge(v, v, G);
    }
}

/// Build an initiator matrix, used as the seed for a larger graph. If the large
/// graph does not need to be build, return a minimal graph with two connected vertices.
/// All graph information on the Kronecker product graph will be stored in the
/// first vertex.
template<typename Graph>
Graph build_initiator_graph (const bool entire_graph_needed) {
    if (entire_graph_needed){
        Graph G{1};
        Utils::add_self_edges(G);
        return G;
    }
    else {
        Graph G{2};
        boost::add_edge(0, 1, G);
        return G;
    }
}
// .............................................................................
// ... Property extraction and analysis functions ..............................
// .............................................................................

/// Collect all parameters that are either target parameters or merely to be
/// calculated from a configuration. Returns a map of keys (parameter names),
/// whether they are to be calculated (true/false + value), and whether they are
/// targets (true/false + value).
/**
  * \param analysis_cfg        The cfg of parameters to be analysed
  * \param target_cfg          The cfg of target parameters

  * \return analysis_targets   A map of parameter names as keys, and a second map
                               as the value, which contains the values for the
                               calculations, as well as the target values
  */
map_type get_analysis_targets(const Config& analysis_cfg,
                              const Config& target_cfg = YAML::Node())
{
    using namespace std;

    // Number of vertices are always calculated for Kronecker products
    map_type analysis_targets
    = {{"num_vertices", entry_type{{"calculate", make_pair<const bool, double>(true, 1)}}}};

    // First, get potential target parameters. Target parameters are always
    // automatically calculated.
    if (target_cfg.size()>0) {
        analysis_targets["error"]
          = entry_type{{"calculate", make_pair<const bool, double>
                       (true, 0)
                      }};
        analysis_targets["largest_comp"]
          = entry_type{{"calculate", make_pair<const bool, size_t>
                       (true, 0)
                      }};

        analysis_targets["num_factors"]
          = entry_type{{"calculate", make_pair<const bool, size_t>
                       (true, 0)
                      }};

        analysis_targets["num_Paretos"]
          = entry_type{{"calculate", make_pair<const bool, size_t>
                       (true, 0)
                      }};

        analysis_targets["num_vertices"]["target"]
          = make_pair<const bool, const double>(true, get_as<double>("value", target_cfg["num_vertices"]));

        for (const auto& param : target_cfg){
            const string param_name = param.first.as<string>();
            if (param_name == "degree_sequence") {

                if (get_as<string>(param_name, target_cfg) == ""){
                    continue;
                }
                analysis_targets[param_name]
                = entry_type{{"calculate", make_pair<bool, vector_pt>(true, vector_pt{{0, 1}})},
                             {"target", make_pair<const bool, const string>
                             (true, get_as<string>(param_name, target_cfg))
                            }};
            }
            else if (param_name == "num_vertices"){
                analysis_targets[param_name]
                = entry_type{{"calculate", make_pair<bool, double>(true, 1)},
                             {"target", make_pair<const bool, const double>
                             (true, get_as<double>("value", target_cfg[param_name]))
                            }};
            }
            else {
                analysis_targets[param_name]
                = entry_type{{"calculate", make_pair<bool, double>(true, 0)},
                             {"target", make_pair<const bool, const double>
                             (true, get_as<double>("value", target_cfg[param_name]))
                            }};
            }
            if (param_name == "clustering") {
                analysis_targets["degree_variance"]
                = entry_type{{"calculate", make_pair<bool, double>(true, 0)},
                             {"target", make_pair<const bool, const double>(false, -1)
                            }};
            }
        }
    }

    // Next, add any potential analysis targets that are not already included in
    // the target config
    for (const auto& param : analysis_cfg){
        const string param_name = param.first.as<string>();

        if (analysis_targets.find(param_name) == analysis_targets.end()) {


            if (get_as<bool>(param_name, analysis_cfg) == false) {
              continue;
            }
            if (param_name == "enabled") {
              continue;
            }

            if (param_name == "degree_sequence") {
                analysis_targets[param_name]
                = entry_type{{"calculate", make_pair<bool, vector_pt>(true, vector_pt{{0, 1}})},
                             {"target", make_pair<const bool, const std::string>
                             (false, "")
                            }};
            }
            else {
                analysis_targets[param_name]
                = entry_type{{"calculate", make_pair<const bool, double>(true, 0)},
                             {"target", make_pair<const bool, const int>(false, -1)}};
            }

            // If the clustering coefficient is included either as a target or
            // analysis parameter, the mean degrees and the degree variance must
            // also be calculated
            if (param_name == "clustering") {
                if (analysis_targets.find("mean_degree") == analysis_targets.end()){
                    analysis_targets["mean_degree"]
                      = entry_type{{"calculate", make_pair<const bool, double>(true, 0)},
                                   {"target", make_pair<const bool, const double>(false, -1)}};
                }
                if (analysis_targets.find("degree_variance") == analysis_targets.end()){
                    analysis_targets["degree_variance"]
                      = entry_type{{"calculate", make_pair<const bool, double>(true, 0)},
                                   {"target", make_pair<const bool, const double>(false, -1)}};
                }
            }
        }
    }

    return analysis_targets;

}

/// Calculates the graph properties with a given map of analysis parameters
/**
  * \tparam Graph              The graph type
  * \tparam Logger             The model logger type
  *
  * \param h                   The graph to be analysed
  * \param analysis_targets    The map of calculated parameters from previous
                               factors
  * \param log                 The model logger
  */
template<typename Graph, typename Logger>
void calculate_properties(const Graph& h,
                          map_type& analysis_targets,
                          const Logger& log)
{
    for (auto& p : analysis_targets) {
        const std::string param = p.first;
        if (analysis_targets[param]["calculate"].first == false){
            continue;
        }

        log->debug("Calculating parameter {} ...", param);

        double k, v = 0;
        std::pair<double, double> degree_statistics;
        if ((analysis_targets.find("degree_variance") != analysis_targets.end())
            or (analysis_targets.find("mean_degree") != analysis_targets.end())
            or (analysis_targets.find("clustering") != analysis_targets.end()))
        {
            degree_statistics = NetworkAnalyser::degree_statistics(h);
        }

        const auto val = analysis_targets[param]["calculate"].second;

        if (param == "clustering") {
            k = std::any_cast<double>(analysis_targets["mean_degree"]["calculate"].second);
            v = std::any_cast<double>(analysis_targets["degree_variance"]["calculate"].second);
            analysis_targets[param]["calculate"].second
            = Kronecker_clustering(std::any_cast<double>(val),
                                   NetworkAnalyser::global_clustering_coeff(h),
                                   k, degree_statistics.first,
                                   v, degree_statistics.second);
        }
        else if (param == "degree_sequence"){
            analysis_targets[param]["calculate"].second
            = Kronecker_degree_sequence(std::any_cast<vector_pt>(val), NetworkAnalyser::degree_sequence(h));
        }
        else if (param == "degree_variance"){
            k = std::any_cast<double>(analysis_targets["mean_degree"]["calculate"].second);
            analysis_targets[param]["calculate"].second
            = Kronecker_degree_variance(k,
                                        degree_statistics.first,
                                        std::any_cast<double>(val),
                                        degree_statistics.second);
        }
        else if (param == "diameter") {
            analysis_targets[param]["calculate"].second
            = std::max(std::any_cast<double>(val), NetworkAnalyser::diameter(h));
        }
        else if (param == "largest_comp") {
            analysis_targets[param]["calculate"].second
            = std::max(std::any_cast<size_t>(val), boost::num_vertices(h));
        }
        else if (param == "mean_degree") {
            analysis_targets[param]["calculate"].second
            = Kronecker_mean_degree(std::any_cast<double>(val), degree_statistics.first);
        }
        else if (param == "num_factors") {
            analysis_targets[param]["calculate"].second = std::any_cast<size_t>(val)+1;
        }
        else if (param == "num_vertices"){
            analysis_targets[param]["calculate"].second
            = Kronecker_num_vertices(std::any_cast<double>(val), boost::num_vertices(h));
        }
    }
}

/// Writes the graph properties into the first vertex of the resulting graph
/**
  * \param K                   The final Kronecker graph
  * \param analysis_targets    The map of calculated parameters from previous factors
  */
template<typename Graph>
void write_properties(Graph& K, map_type& analysis_targets){
    for (auto& p : analysis_targets) {
        const std::string param = p.first;
        if (analysis_targets[param]["calculate"].first == false){
            continue;
        }
        auto val = analysis_targets[param]["calculate"].second;

        if (param == "clustering") {
            K[0].state.clustering_global = std::any_cast<double>(val);
        }
        else if (param == "degree_sequence") {
            K[0].state.degree_sequence = std::any_cast<vector_pt>(val);
        }
        else if (param == "degree_variance") {
            K[0].state.degree_variance = std::any_cast<double>(val);
        }
        else if (param == "diameter") {
            K[0].state.diameter = std::any_cast<double>(val);
        }
        else if (param == "error") {
            K[0].state.error = std::any_cast<double>(val);
        }
        else if (param == "largest_comp") {
            K[0].state.largest_comp = std::any_cast<size_t>(val);
        }
        else if (param == "mean_degree") {
            K[0].state.mean_degree = std::any_cast<double>(val);
        }
        else if (param == "num_factors") {
            K[0].state.num_factors = std::any_cast<size_t>(val);
        }
        else if (param == "num_Paretos") {
            K[0].state.num_Paretos = std::any_cast<size_t>(val);
        }
        else if (param == "num_vertices") {
            K[0].state.num_vertices = std::any_cast<double>(val);
        }
    }
}

// Extracts just the target names and numerical values from a map of analysis
// and target properties
std::map<std::string, double> extract_targets(map_type& analysis_and_targets)
{
    std::map<std::string, double> targets;
    for (auto& p : analysis_and_targets) {
        const std::string param = p.first;
        if (analysis_and_targets[param]["target"].first == false) {
            continue;
        }
        if (param == "degree_sequence") {
            continue;
        }

        targets[param]
        = std::any_cast<double>(analysis_and_targets[param]["target"].second);

    }

    return targets;
}

/// Return 'true' if a parameter is a target/analysis parameter
bool is_param(map_type& map, std::string param, std::string target = "target") {

    if (map.find(param) != map.end()){
        if (map[param][target].first) {
            return true;
        }
    }
    return false;
}
/// Overloads
/// Return 'true' if a parameter is a target and is equal to a given value
bool is_target(map_type& map, std::string param, std::string val){
    if (map.find(param) != map.end()){
        if (map[param]["target"].first) {
            if (std::any_cast<std::string>(map[param]["target"].second) == val) {
                return true;
            }
        }
    }
    return false;
}

bool is_target(map_type& map, std::string param, std::size_t val){
    if (map.find(param) != map.end()){
        if (map[param]["target"].first) {
            if (std::any_cast<size_t>(map[param]["target"].second) == val) {
                return true;
            }
        }
    }
    return false;
}

bool is_target(map_type& map, std::string param, double val){
    if (map.find(param) != map.end()){
        if (map[param]["target"].first) {
            if (std::any_cast<double>(map[param]["target"].second) == val) {
                return true;
            }
        }
    }
    return false;
}

// Find an estimated mu value for a given clustering coefficient
// TO DO: This must be generalised to arbitrary target parameters
double mu_inv(const double N, const double k, const double c_t)
{
    double b = 1.8;
    double m = 1.0*std::round(-0.5*sqrt(4.0*(pow(N, 2)) - 4.0*N*(k+1) + 1) + N - 0.5);
    double c_0 = 2./3.-1/(3*m);
    double c_1 = pow((m-1), (2./3))/(3.3)*pow(log(N), 2)/N;
    double a = (c_0 - c_1)*exp(b);

    return std::max(0., std::min(1., sqrt(-1.0/b*log((c_t-c_1)/a))-1));
}

// .............................................................................
// ... Grid search utility functions ...........................................
// .............................................................................

// Returns pairs of positive integer factors, including (1, N) if specified
// Throws an error if N is not strictly positive
factors get_factors (const size_t N, const bool include_trivial=false) {

    if (N <= 0) {
        throw std::invalid_argument("N must be positive!");
    }

    factors res {{1, N}};
    if (not include_trivial){
        res.pop_back();
    }
    const size_t limit = int(sqrt(N)+1);
    res.reserve(limit);

    for (size_t i = 2; i < limit; ++i) {
        if (not (N%i)) {
            const factor t = {i, static_cast<size_t>(N/i)};
            res.emplace_back(t);
        }
    }

    return res;
}

// Returns pairs of positive integers such that (i+1)(j+1)=(k+1),
// including (0, k) if specified
factors get_k_factors (const size_t k, const bool include_trivial=false) {

    factors res {{0, k}};
    if (not include_trivial){
        res.pop_back();
    }
    res.reserve(k);

    for (size_t i = 1; i < k; ++i) {
        if (not ((k+1)%(i+1))) {
            const factor t = {i, static_cast<size_t>((k+1)/(i+1)-1)};
            res.emplace_back(t);
        }
    }

    return res;
}

// Factorizes a list of d-tuple factors into a list of d+1-tuple factors
factors factorize(const factors& target, factors (*factor_func)(size_t, bool)) {

    factors res = {};

    for (const auto& fac: target) {
        for (size_t i=0; i < fac.size(); ++i) {
            const auto x = factor_func(fac[i], false);
            for (const auto& x_fac: x) {
                factor q = factor(fac.begin(), fac.begin()+i);
                q.insert(q.end(), x_fac.begin(), x_fac.end());
                q.insert(q.end(), fac.begin()+(i+1), fac.end());
                sort(q.begin(), q.end());
                bool already_in = false;
                for (const auto& c : res) {
                    if (c == q) {
                        already_in = true;
                        break;
                    }
                }
                if (not already_in) {
                    res.emplace_back(q);
                }
            }
        }
    }

    return res;
}

// Returns d-tuples of positive integer factors of N
/*
 * \param val                   The natural number to factorise
 * \param factor_func           The factorisation function
 * \param d                     The factor size
 *
 * \returns res                 A list of tuples
 * \throws invalid_argument     If d_min > d_max
 */
factors get_d_factors (const size_t val,
                       factors (*factor_func)(size_t, bool),
                       const size_t d)
{
    if (d < 1) {
        throw std::invalid_argument("d must be greater than 1!");
    }
    if (d == 1) {
        return {{val}};
    }

    factors res = factor_func(val, false);
    for (size_t p = 2; p < d; ++p) {
        res = factorize(res, factor_func);
    }

    return res;

}

// Get a grid of d-factors (including reducible factors) for a given error range
/*
 * \param target        The grid center
 * \param factor_func   The factorisation function
 * \param err           The maximum grid boundary, given as a relative error of
                        the grid center
 * \param d             The max factorisation length
 */
factors get_grid (const size_t grid_center,
                  factors (*factor_func)(size_t, bool),
                  const double err,
                  const size_t d_min,
                  const size_t d_max)
{
    factors res = {};
    const auto delta = int(grid_center*err);
    for (size_t i = grid_center - delta; i <= grid_center + delta; ++i) {
        for (size_t p = d_min; p <= d_max; ++p) {
            const auto t = get_d_factors(i, factor_func, p);
            if (t.size() == 0){
              continue;
            }
            res.insert(res.end(), t.begin(), t.end());
        }
    }

    return res;
}

// Return the N-grid (using the standard factorisation)(single dimension)
factors get_N_grid (const size_t grid_center,
                    const double err,
                    const size_t d)
{
    return get_grid(grid_center, get_factors, err, d, d);
}

// Return the N-grid (using the standard factorisation)(dimension span)
factors get_N_grid (const size_t grid_center,
                    const double err,
                    const size_t d_min,
                    const size_t d_max)
{
    return get_grid(grid_center, get_factors, err, d_min, d_max);
}

// Return the k-grid
factors get_k_grid (const size_t grid_center,
                    const double err,
                    const size_t d)
{
    return get_grid(grid_center, get_k_factors, err, d, d);
}

// Return the k-grid (dimension span)
factors get_k_grid (const size_t grid_center,
                    const double err,
                    const size_t d_min,
                    const size_t d_max)
{
    return get_grid(grid_center, get_k_factors, err, d_min, d_max);
}

// Return the mu-grid
std::vector<double> get_mu_grid(const double mu,
                                const double tolerance,
                                const size_t grid_size)
{
    const double lower_bound = std::max(0., mu*(1-tolerance));
    const double upper_bound = std::min(1., mu*(1+tolerance));
    const double step_size = (upper_bound-lower_bound)/grid_size;
    std::vector<double> grid;
    grid.reserve(grid_size);
    for (size_t i = 0; i < grid_size; ++i) {
        grid.emplace_back(lower_bound + i*(step_size));
    }

    return grid;
}

} // namespace KronGen::Utils

#endif // UTOPIA_MODELS_KRONGEN_UTILS
