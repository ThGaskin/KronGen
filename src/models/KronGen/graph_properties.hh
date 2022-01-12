#ifndef UTOPIA_MODELS_KRONGEN_GRAPHPROPERTIES
#define UTOPIA_MODELS_KRONGEN_GRAPHPROPERTIES

#include "graph_types.hh"

namespace Utopia::Models::KronGen::GraphProperties {

using namespace Utopia::Models::KronGen::GraphTypes;

// ... Estimators for clustering coefficients  .................................

/// Returns the clustering coefficient of an Erdos-Renyi random graph
double clustering_ER (const double& N, const double& k)
{

    // undefined, but resulting Kronecker c does not depend on this
    if ((N == 2 and k == 1) or (N == 1)) {
        return 0;
    }

    return (k/(N-1));
}

/// Returns the clustering coefficient of regular graph
double clustering_regular(const size_t& N, const size_t& k)
{
    if (k < 2) {
        return 0;
    }
    if (N <= k+1) {
        return 1;
    }
    else {
        double res = 3./4.*(pow(k, 2)-2*k)/(pow(k, 2)-k);
        // Additional triangles through looping over network
        if (3.0*k/2 >= 1.0*N) {
            auto a = ((3*k/2)%N+1)*((3*k/2)%N+2);
            double b = (static_cast<double>(a))/(k*k-1);
            res += b;
        }
        return res;
    }
}

/// Returns an estimate of the clusering coefficient of various graph types
/**
 * \param N     The number of vertices
 * \param k     The mean degree
 * \param t     The graph type
 *
 * \return c    The estimated clustering coefficient
 */
double clustering_estimation(const double& N, const double& k, const GraphType& t)
{
    if (t == GraphType::Chain) {
        return 0;
    }
    else if (t == GraphType::Complete) {
        return 1;
    }
    else if (t == GraphType::Regular){
        return clustering_regular(N, k);
    }
    else {
        return clustering_ER(N, k);
    }
}

// ... Estimators for diameters ................................................

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
// To do: factor of 1.7?
double diameter_estimation_ER (const double& N, const double& k)
{
    if (k == 0) {
        return -1;
    }
    else if (k == 1) {
        return (N-1);
    }
    else if (k == N-1) {
        return 1;
    }
    else if (N == 2) {
        return 1;
    }
    else {
        return (1.7 * log(N)/log(k));
    }
}

/// Returns an estimation of the diameter of a regular network, given the number of
/// vertices and the mean degree
// To do: need better estimate here
double diameter_estimation_regular (const double& N, const double& k)
{
    return ceil(floor(N/2) / (k/2));
}

/// Returns an estimate of the diameter of various graph types
/**
 * \param N     The number of vertices
 * \param k     The mean degree
 * \param t     The graph type
 *
 * \return d    The estimated diameter
 */
double diameter_estimation(const double& N, const double& k, const GraphType& t)
{
    if (t == GraphType::Chain) {
        return (N-1);
    }
    else if (t == GraphType::Regular){
        return diameter_estimation_regular(N, k);
    }
    else {
        return diameter_estimation_ER(N, k);
    }
}

// ... Mean degrees ............................................................

/// Returns the mean degree of a 1-d chain of length N
double mean_degree_chain (const double& N)
{
    return (2*(N-1)/N);
}

/// Returns the mean degree various graph types
/**
 * \param N     The number of vertices
 * \param k     The input mean degree
 * \param t     The graph type
 *
 * \return d    The actual mean degree
 */
double mean_degree (const double& N, const double& k, const GraphType& t)
{
    if (t == GraphType::Chain) {
        return mean_degree_chain(N);
    }
    else if (t == GraphType::Complete){
        return (N-1);
    }
    else {
        return k;
    }
}

// ... Degree variance .........................................................

/// Returns the degree variance of a 1-d chain of length N
double degree_variance_chain (const double& N)
{
    const double k = mean_degree_chain(N);

    return (1/N * (2*pow((1-k), 2) + (N-2)*pow((2-k), 2)));
}

/// Returns the degree variance of an ErdosRenyi graph
double degree_variance_ER (const double N, const double k)
{
    if (N == 1) {
        return 0;
    }
    else {
        return k*(1-k/(N-1));
    }
}

/// Returns the degree variance of various graph types
/**
 * \param N     The number of vertices
 * \param k     The mean degree
 * \param t     The graph type
 *
 * \return d    The degree variance
 */
// To do: other graph types, e.g. scale-free
double degree_variance(const double& N, const double& k, const GraphType t)
{
    if (t == GraphType::Chain){
        return degree_variance_chain(N);
    }
    else if (t == GraphType::Complete) {
        return 0;
    }
    else if (t == GraphType::ErdosRenyi) {
        return degree_variance_ER(N, k);
    }
    else if (t == GraphType::Regular) {
        return 0;
    }
    else {
        return -1;
    }
}

} // namespace KronGen::GraphProps

#endif // UTOPIA_MODELS_KRONGEN_GRAPHPROPERTIES
