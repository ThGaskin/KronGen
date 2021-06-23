#ifndef UTOPIA_MODELS_KRONGEN_UTILS
#define UTOPIA_MODELS_KRONGEN_UTILS

#include <random>
#include <algorithm>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include "utopia/core/graph/iterator.hh"

#include "../NetworkAnalyser/graph_metrics.hh"

namespace Utopia::Models::KronGen::Utils {

using namespace std;
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
    using namespace Utopia;

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

/// Caclulate the required mean degree of a Kronecker factor, given the target
/// and one other factor
double Kronecker_mean_degree_factor (const double m, const double g)
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
  * \param h        The degree distribution variance of H
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
    using vertices_size_type = typename boost::graph_traits<Graph>::vertices_size_type;

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
        const auto starting_point = fourSweep<vertices_size_type>(h);
        const double d = iFUB(starting_point.first, starting_point.second, 0, h);
        diam = max(diam, d);
    }
}


// ... Utility functions .......................................................

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
// a desired product
vector_pt N_factors (const size_t N)
{
    vector<bool> candidates(round(static_cast<double>(N)/2.+2), true);
    vector_pt res;
    for (int i = 2; i < round(static_cast<double>(N)/2+1); ++i) {
        if (not (N % i) && (candidates[i]==true)) {
            res.push_back({i, N/i});
            candidates[N/i] = false;
        }
    }
    return res;
}

/// Return a list of possible vertex number factor pairs producing a desired
/// product, or the closest possible factor pairs
vector_pt closest_N_factors (const size_t N)
{
    auto res = N_factors(N);
    if (res.empty()) {
        res = closest_N_factors(N+1);
        if (N>4) {
            auto res2 = closest_N_factors(N-1);
            res.insert(res.end(), res2.begin(), res2.end());
        }
    }

    return res;
}

/// Return a list of all possible mean degree factors producing a desired product
vector_pt mean_deg_factors (const size_t m)
{
    vector<bool> candidates(round(static_cast<double>(m)/2.+2), true);
    vector_pt res;
    for (int i = 1; i < round(static_cast<double>(m)/2); ++i) {
        if (not ((m+1) % (i+1)) && (candidates[i]==true)) {
            auto k = ((m+1)/(i+1)-1);
            res.push_back({i, k});
            candidates[k] = false;
        }
    }

    return res;
}

/// Return a list of possible mean degree number factor pairs producing a desired
/// product, or the closest possible factor pairs
vector_pt closest_mean_deg_factors (const size_t m)
{
    auto res = mean_deg_factors(m);
    if (res.empty()) {
        res = closest_mean_deg_factors(m+1);
        if (m>3) {
            auto res2 = closest_mean_deg_factors(m-1);
            res.insert(res.end(), res2.begin(), res2.end());
        }
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

/// Returns an estimation of the diameter of a network, given the number of
/// vertices and the mean degree
/**
  * \param N     The number of vertices
  * \param m     The mean degree
  *
  * \return d    The estimated diameter
*/
// To do: test me
// To do: extreme cases
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
        return (1.6 * log(N)/log(m));
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

}
#endif
