#ifndef UTOPIA_MODELS_KRONGEN_UTILS
#define UTOPIA_MODELS_KRONGEN_UTILS

#include <algorithm>
#include <random>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include "utopia/core/graph/iterator.hh"

#include "../NetworkAnalyser/graph_metrics.hh"

namespace Utopia::Models::KronGen::Utils {

using namespace std;
using factor = typename std::vector<std::size_t>;
using factors = typename std::vector<std::vector<std::size_t>>;
using pair_pt = typename std::pair<std::pair<size_t, double>, std::pair<size_t, double>>;
using vector_pt = typename std::vector<pair<size_t, size_t>>;

// ... Kronecker properties ....................................................
/// Kronecker product of graphs. Graphs must have a self-loop on every node
/**
  * \param G    The first factor
  * \param H    The second factor
  *
  * \return K   The Kronecker product
*/
template<typename Graph>
Graph Kronecker_product(Graph& G, Graph& H) {

    using namespace boost;

    const size_t M = num_vertices(H);
    Graph K{num_vertices(G) * M};

    for(const auto g : Utopia::range<IterateOver::edges>(G)) {
        for(const auto h : Utopia::range<IterateOver::edges>(H)) {
            const auto s_1 = source(g, G);
            const auto s_2 = source(h, H);
            const auto t_1 = target(g, G);
            const auto t_2 = target(h, H);

            auto i = s_1*M + s_2;
            auto j = t_1*M + t_2;
            if ((not edge(i, j, K).second) && (not edge(j, i, K).second)){
                add_edge(i, j, K);
            }

            i = s_1*M + t_2;
            j = t_1*M + s_2;
            if ((not edge(i, j, K).second) && (not edge(j, i, K).second)){
                add_edge(i, j, K);
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

/// Calculates the required mean degree of a graph in order to generate a given
/// target clustering coefficient. The graph can either be complete or have
/// zero clustering.
/**
  * \param is_zero  Whether the graph in question has zero clustering or not
  * \param c_H      The clustering coefficient of the other Kronecker factor H
  * \param c        The target clustering coefficient
  * \param h        The mean degree of graph H
  * \param v        The degree distribution variance of H
*/
double get_mean_deg_c(const double c_H,
                      const double c,
                      const double h,
                      const double v)
{
    double g;

    // one Kronecker factor has clustering coefficient zero
    if (c_H >= c) {
        g = 3*pow(h, 2)*c_H-2*pow(h, 2)*c-3*h*c_H-h*c+6*h+3*v*c_H-2*c*v+c;
        g = pow(g, 2);
        g += -4*(-pow(h, 2)*c-2*h*c-c*v-c)*(pow(h, 2)*c_H-pow(h, 2)*c-h*c_H+h*c+c_H*v-c*v);
        g = sqrt(g);
        g += 3*pow(h, 2)*c_H -2*pow(h, 2)*c-3*h*c_H-h*c+6*h+3*c_H*v-2*c*v+c;
        g /= 2*c*(pow(h, 2)+2*h+v+1);
    }

    // one Kronecker factor is complete
    else {
        g = 8*pow(h, 2)*c_H+pow(h, 2)*c-9*pow(h, 2)-8*h*c_H+2*h*c+6*h+8*c_H*v-8*c*v+c-1;
        g *= c-1;
        g = sqrt(g);
        g += -2*pow(h, 2)*c_H+2*pow(h, 2)*c+2*h*c_H+h*c-3*h-2*c_H*v+2*c*v-c+1;
        g /= (2*(pow(h, 2)*c_H-pow(h, 2)*c-h*c_H-2*h*c+3*h+c_H*v-c*v-c+1));
    }

    // always select positive root
    return (fabs(g));

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
    using namespace Utopia::Models::NetworkAnalyser;

    // Calculate the clustering coefficient
    if (calculate_c) {
        const double c_temp = global_clustering_coeff(h);
        const auto deg_stats = degree_statistics(h);
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
        diam = max(diam, diameter(h));
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

// Return the N-grid (using standard factorisation)
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

// Error norm
double err_func(double x, double y) {
    return std::abs(1-x/y); // L1 norm
    //return pow((1-x/y), 2); //L2 norm
}

// Objective function for N
double N_err (double N_t, factor N, factor k) {
    double N_res = 1.0;
    for (const auto& n : N) {
        N_res *= n;
    }

    return err_func(N_res, N_t);
}

// Objective function for k
double k_err (double k_t, factor N, factor k) {
    double k_res = 1.0;
    for (const auto& kk : k) {
        k_res *= (kk+1);
    }

    return err_func(k_res, k_t);
}

// ... General utility functions ...............................................

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

/// Return a list of possible vertex number factor pairs producing
// a desired product. Does not include (1, N) or (2, N/2)
vector_pt N_factors (const size_t N)
{
    vector<bool> candidates(round(static_cast<double>(N)/2.+2), true);
    candidates[2]=false;
    vector_pt res;
    for (int i = 2; i < round(static_cast<double>(N)/2+1); ++i) {
        if (not (N % i) && (candidates[i]==true)) {
            if (candidates[N/i] == true) {
                res.push_back({i, N/i});
                candidates[N/i] = false;
            }
        }
    }
    return res;
}

/// Return a list of possible vertex number factor pairs producing a desired
/// product, or the closest possible factor pairs
vector_pt closest_N_factors (const size_t N)
{
    auto res = N_factors(N);
    size_t i = 1;
    while (res.empty()) {
        auto res2 = N_factors(N+i);
        res.insert(res.end(), res2.begin(), res2.end());
        if (N > i+2){
            res2=N_factors(N-i);
            res.insert(res.end(), res2.begin(), res2.end());
        }
        ++i;
    }
    return res;
}

/// Return a list of all possible mean degree factors producing a desired product
vector_pt mean_deg_factors (const size_t m)
{
    vector<bool> candidates(round(static_cast<double>(m)/2.+2), true);
    vector_pt res;
    candidates[1]=false;
    for (int i = 2; i < round(static_cast<double>(m)/2); ++i) {
        if (not ((m+1) % (i+1)) && (candidates[i]==true)) {
            auto k = ((m+1)/(i+1)-1);
            if (candidates[k] == true) {
                res.push_back({i, k});
                candidates[k] = false;
            }
        }
    }

    return res;
}

/// Return a list of possible mean degree number factor pairs producing a desired
/// product, or the closest possible factor pairs
vector_pt closest_mean_deg_factors (const size_t m)
{
    auto res = mean_deg_factors(m);
    size_t i = 1;
    while (res.empty()) {
        auto res2 = mean_deg_factors(m+i);
        res.insert(res.end(), res2.begin(), res2.end());
        if (m > i+3){
            res2=mean_deg_factors(m-i);
            res.insert(res.end(), res2.begin(), res2.end());
        }
        ++i;
    }

    return res;
}

/// Return two factors of a pair consisting of {num_vertices, mean_degree}
template<typename RNGType>
pair_pt get_factors_N_m (const size_t N,
                         const size_t m,
                         RNGType& rng)
{
    // Get the Kronecker factors of N, m and shuffle
    auto N_factors = Utils::closest_N_factors(N);
    auto m_factors = Utils::closest_mean_deg_factors(m);
    shuffle(N_factors.begin(), N_factors.end(), rng);
    shuffle(m_factors.begin(), m_factors.end(), rng);

    for (const auto& N_i : N_factors) {
        for (const auto& m_i : m_factors) {
            if ((N_i.first > m_i.first) and (N_i.second > m_i.second)) {
                return {{N_i.first, m_i.first}, {N_i.second, m_i.second}};
            }
            else if ((N_i.second > m_i.first) and (N_i.first > m_i.second)) {
                return {{N_i.second, m_i.first}, {N_i.first, m_i.second}};
            }
        }
    }

    return {{0, 0},{0, 0}};
}

/// Returns an estimation of the diameter of a random network, given the number of
/// vertices and the mean degree
/**
  * \param N     The number of vertices
  * \param m     The mean degree
  *
  * \return d    The estimated diameter
*/
// To do: test me
// To do: extreme cases
// To do: factor of 1.6?
double diameter_estimation (const double N, const double m) {
    if (m == 0) {
        return -1;
    }
    else if (m == 1) {
        return (N-1);
    }
    else if (N == 2) {
        return 1;
    }
    else {
        return (1.7 * log(N)/log(m));
    }
}

/// Returns the mean degree of a 1-d chain of length N
/**
  * \param N    The length of the chain
  *
  * \return m   The mean degree of the graph
*/
double mean_degree_chain_graph (const double N) {

    return (2*(N-1)/N);
}

/// Returns the clustering coefficient of Erdos-Renyi random graph
double ER_clustering (const double N, const double k) {

    if (N == 1) {return 0;} // undefined, but resulting Kronecker c does not depend on this
    return (k/(N-1));
}

double ER_variance (const double N, const double k) {
    if (N == 1) {return 0;}
    return k*(1-k/(N-1));
}

/// Returns the clustering coefficient of regular graph
/**
  * \param N    The number of vertices
  * \param m    The mean degree
  *
  * \return c   The clustering coefficient
*/
// To do: linear interpolation not quite accurate
double regular_graph_clustering(const double& N, const double& m) {
    if (N < 3 or m < 2) {
        return 0;
    }
    double T = (m/2 * (m-1.))-m/4*(m/2+1);
    const double T_2 = (m-2*N/3)*(m/4*(m/2+1))/(N/3-1.);
    if ( 2*N/3 < m) {
        T += T_2;
    }

    return (T /(m/2*(m-1)));
}

}
#endif
