#ifndef UTOPIA_MODELS_KRONGEN_OBJECTIVE_FUNCS
#define UTOPIA_MODELS_KRONGEN_OBJECTIVE_FUNCS
#define UNUSED(expr) do { (void)(expr); } while (0)

#include "graph_properties.hh"
#include "graph_types.hh"
#include "utils.hh"

namespace Utopia::Models::KronGen::ObjectiveFuncs {

using namespace Utopia::Models::KronGen;
using namespace Utopia::Models::KronGen::GraphTypes;

using factor = typename std::vector<std::size_t>;

// Error norm
// To do: initialise this from the front end
double err_func(double x, double y) {
    return std::abs(1-x/y); // L1 norm
    //return pow((1-x/y), 2); //L2 norm
}

/// Clustering objective function for various graphs types
/**
  * \param c_t     The target clustering coefficient
  * \param N       Vector of graph vertex counts
  * \param k       Vector of graph mean degrees
  * \param t       Vector of graph types
  *
  * \return err    The error in the specified norm
  *
  * @throws        Invalid argument if number of factors do not match
  */
// Missing: clustering for scale-free graphs etc.
double clustering_obj_func (const double& c_t,
                            const factor& N,
                            const factor& k,
                            const std::vector<GraphType>& t)
{
    if ((N.size() != k.size()) or (k.size() != t.size())) {
        throw std::invalid_argument("Number of factors do not match!");
    }

    // Calculate graph degree variances and clustering coefficients
    double c_current=0., c_previous, k_current, k_previous, v_current, v_previous;

    for (size_t i = 0; i < N.size(); ++i){
        k_current = k[i];
        c_current = GraphProperties::clustering_estimation(N[i], k[i], t[i]);
        v_current = GraphProperties::degree_variance(N[i], k[i], t[i]);
        if (i > 0) {
            c_current = Utils::Kronecker_clustering(c_current, c_previous,
                                                    k_current, k_previous,
                                                    v_current, v_previous);
            v_current = Utils::Kronecker_degree_variance(k_current, k_previous,
                                                         v_current, v_previous);
            k_current = Utils::Kronecker_mean_degree(k_current, k_previous);
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
  * \param N       Vector of graph vertex counts
  * \param k       Vector of graph mean degrees
  * \param t       Vector of graph types
  *
  * \return err    The error in the specified norm
  *
  * @throws        Invalid argument if number of factors do not match
  */
// Missing: diameter for scale-free graphs etc.
double diameter_obj_func (const double& d_t,
                          const factor& N,
                          const factor& k,
                          const std::vector<GraphType>& t)
{
    if ((N.size() != k.size()) or (k.size() != t.size())) {
        throw std::invalid_argument("Number of factors do not match!");
    }

    double d = 0;
    for (size_t i = 0; i < N.size(); ++i){
        const size_t d_est = GraphProperties::diameter_estimation(N[i], k[i], t[i]);
        if (d_est > d) {
            d = d_est;
        }
    }

    return err_func(d, d_t);

}

/// Number of vertices objective function, independent of graph type
/**
  * \param d_t     The target vertex count
  * \param N       Vector of graph vertex counts
  * \param k       Vector of graph mean degrees
  * \param t       Vector of graph types
  *
  * \return err    The error in the specified norm
  */
 double N_obj_func (const double& N_t,
                    const factor& N,
                    const factor& k,
                    const std::vector<GraphType>& t)
{
    UNUSED(k);
    UNUSED(t);
    double N_res = 1;
    for (const auto& n : N) {
        N_res *= n;
    }

    return err_func(N_res, N_t);
}

/// Mean degree objective function for various graph types
/*
 * \param d_t     The target mean degree
 * \param N       Vector of graph vertex counts
 * \param k       Vector of graph mean degrees
 * \param t       Vector of graph types
 *
 * \return err    The error in the specified norm
 */
double k_obj_func (const double& k_t,
                   const factor& N,
                   const factor& k,
                   const std::vector<GraphType>& t)
{
    double k_res = 1.0;
    for (size_t i = 0; i < k.size(); ++i) {
        k_res *= (GraphProperties::mean_degree(N[i], k[i], t[i])+1);
    }

    return err_func(k_res, k_t+1);
}

} // namespace KronGen::ObjectiveFuncs

#endif // UTOPIA_MODELS_KRONGEN_OBJECTIVE_FUNCS
