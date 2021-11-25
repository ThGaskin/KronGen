#ifndef UTOPIA_MODELS_KRONGEN_UTILS
#define UTOPIA_MODELS_KRONGEN_UTILS

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include "utopia/core/graph/iterator.hh"

#include "../NetworkAnalyser/graph_metrics.hh"

namespace Utopia::Models::KronGen::Utils {

using factor = typename std::vector<size_t>;
using factors = typename std::vector<std::vector<size_t>>;
using pair_pt = typename std::pair<std::pair<size_t, double>, std::pair<size_t, double>>;
using vector_pt = typename std::vector<std::pair<size_t, size_t>>;

// ... Kronecker properties ....................................................
/// Kronecker product of graphs. Graphs must have a self-loop on every node
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
size_t Kronecker_num_vertices (const size_t N, const size_t M)
{
    return (N*M);
}

/// Calculate the mean degree of a Kronecker product of two graphs G, H
double Kronecker_mean_degree (const double g, const double h)
{
    return ((g+1)*(h+1)-1);
}

/// Calculate the degree distribution variance of a Kronecker product of two
/// graphs G, H with mean degrees g, h, and degree distribution variances
/// v_g, v_h
double Kronecker_degree_variance (const double g,
                                  const double h,
                                  const double v_g,
                                  const double v_h)
{
    return (v_g*v_h + pow(1+g, 2)*v_h + pow(1+h, 2)*v_g);
}

/// Calculate the degree distribution variance of an Erdos-Renyi random graph
/// with N vertices and mean degree k
double ER_degree_variance (const size_t N, const double k) {

    return k*(1-k/(N-1));
}

/// Caclulate the required mean degree of a Kronecker factor, given the target
/// and one other factor
double Kronecker_mean_degree_inv (const double m, const double g)
{
    return ((m+1)/(g+1)-1);
}


/// Calculate the clustering coefficient of a Kronecker product of two graphs G, H
/** This product is symmetric in G, H
 * \param c_G        The clustering coefficient of the first graph
 * \param c_H        The clustering coefficient of the second graph
 * \param g          The mean degree of the first graph
 * \param h          The mean degree of the second graph
 * \param var_g      The degree distribution variance of the first graph
 * \param var_h      The degree distribution variance of the second graph
 *
 * \return c         The clustering coefficient of the Kronecker product
 */
double Kronecker_clustering (const double c_G,
                             const double c_H,
                             const double g,
                             const double h,
                             const double v_g,
                             const double v_h)
{

    const double k = Kronecker_mean_degree(g, h);
    const double var = Kronecker_degree_variance(g, h, v_g, v_h);
    double c = c_G*c_H*(v_g+g*(g-1))*(v_h+h*(h-1));
    c += 3*g*c_H*(v_h+h*(h-1))+3*h*c_G*(v_g+g*(g-1));
    c += c_G*(v_g+g*(g-1))+c_H*(v_h+h*(h-1));
    c += 6*g*h;
    c /= (var + k*(k-1));

    return c;
}

/// Add self_edges to a graph
template<typename Graph>
void add_self_edges (Graph& g) {
    for (const auto& v : range<IterateOver::vertices>(g)) {
        boost::add_edge(v, v, g);
    }
}

/// Remove self_edges from a graph
template<typename Graph>
void remove_self_edges (Graph& g) {
    for (const auto& v : range<IterateOver::vertices>(g)) {
        boost::remove_edge(v, v, g);
    }
}

/// Calculates graph properties (clustering coefficient and diameter) on the go
/// as a Kronecker graph is being generated
template<typename Graph>
void calculate_properties(Graph& h,
                          bool& first_run,
                          const bool& calculate_c,
                          const bool& calculate_diam,
                          double& c_global,
                          double& diam,
                          double& mean_deg,
                          double& variance)
{
    // Calculate the clustering coefficient
    if (calculate_c) {
        const double c_temp = NetworkAnalyser::global_clustering_coeff(h);
        const auto deg_stats = NetworkAnalyser::degree_statistics(h);
        if (first_run) {
            c_global = c_temp;
            mean_deg = deg_stats.first;
            variance = deg_stats.second;
            first_run = false;
        }
        else {
            c_global = Kronecker_clustering(c_global,
                                            c_temp,
                                            mean_deg,
                                            deg_stats.first,
                                            variance,
                                            deg_stats.second);
            variance = Kronecker_degree_variance(mean_deg,
                                                 deg_stats.first,
                                                 variance,
                                                 deg_stats.second);
            mean_deg = Kronecker_mean_degree(mean_deg, deg_stats.first);
        }
    }

    // Calculate the diameter
    if (calculate_diam) {
        diam = std::max(diam, NetworkAnalyser::diameter(h));
    }
}

// ... Grid search utility functions ...........................................
// Returns pairs of positive integer factors, including (1, N) if specified
factors get_factors (const size_t N, const bool include_trivial=true) {

    factors res {{1, N}};
    if (not include_trivial){
        res.pop_back();
    }
    const size_t limit = int(sqrt(N));
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
factors get_k_factors (size_t k, const bool include_trivial=true) {

    factors res {{0, k-1}};
    if (not include_trivial){
        res.pop_back();
    }
    const size_t limit = int(sqrt(k+1));
    res.reserve(limit);

    for (size_t i = 1; i < limit; ++i) {
        if (not ((k+1)%(i+1))) {
            const factor t = {i, static_cast<size_t>((k+1)/(i+1)-1)};
            res.emplace_back(t);
        }
    }

    return res;
}

// Returns d-tuples of positive integer factors of N, including reducible
// tuples (1, 1, 1, ..., q1, ..., qd) if specified
/*
 * \param val           The natural number to factorise
 * \param factor_func   The factorisation function
 * \param d             The number of factors
 * \param include_ones  Whether or not to include reducible tuples (this
                        effectively means including all factors up to d)
 * \returns res         A list of tuples
 */
factors get_d_factors (const size_t val,
                       factors (*factor_func)(size_t, bool),
                       const size_t d_min,
                       const size_t d_max,
                       const bool include_ones)
{
    if (d_max < 2) {
        return {{}};
    }

    auto res = factor_func(val, true);

    if (res.size() == 1 or d_max == 2){
        return res;
    }


    auto current_to_do = res;
    factors current_done = {{}};
    current_to_do.erase(current_to_do.begin());
    for (size_t p = 3; p <= d_max; ++p){
        current_done = {};

        for (const auto& fac: current_to_do) {
            for (size_t i=0; i<fac.size(); ++i) {
                const auto x = factor_func(fac[i], false);
                for (const auto& xx: x) {
                    factor q = factor(fac.begin(), fac.begin()+i); // i-1 ?
                    q.insert(q.end(), xx.begin(), xx.end());
                    q.insert(q.end(), fac.begin()+(i+1), fac.end());
                    sort(q.begin(), q.end());
                    bool already_in = false;
                    for (const auto& c : current_done) {
                        if (c == q) {
                            already_in = true;
                            break;
                        }
                    }
                    if (not already_in) {
                        current_done.emplace_back(q);
                    }
                }
            }
        }
        if (include_ones and p > d_min) {
            res.insert(res.end(), current_done.begin(), current_done.end());
        }
        else {
            res = current_done;
        }
        current_to_do = current_done;
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
    auto delta = int(grid_center*err);
    for (size_t i = grid_center - delta; i <= grid_center + delta; ++i) {
        const auto t = get_d_factors(i, factor_func, d_min, d_max, true);
        res.insert(res.end(), t.begin(), t.end());
    }

    return res;
}

// Return the N-grid (using the standard factorisation)
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
                    const size_t d_min,
                    const size_t d_max)
{
    return get_grid(grid_center, get_k_factors, err, d_min, d_max);
}

} // namespace KronGen::Utils

#endif // UTOPIA_MODELS_KRONGEN_UTILS
