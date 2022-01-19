#ifndef UTOPIA_MODELS_KRONGEN_GRAPHDEF
#define UTOPIA_MODELS_KRONGEN_GRAPHDEF

// third-party library includes
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

namespace Utopia::Models::KronGen::GraphDefinition{

// Defines the graph type used in the KronGen model

// Vertex
struct VertexState
{
    // Betweenness centrality
    double betweenness;

    // Closeness centrality
    double closeness;

    // The local clustering coefficient
    double clustering_local;

    // Core number
    size_t core_number = 0;

    // Degree
    size_t degree;

    // Expected distance to a random other vertex
    double distance_avg;

    // Harmonic average of distances
    double distance_harmonic;

    // Maximum distance to a random other vertex
    double distance_max;

    // Fraction of outgoing links for which the mutual link exists as well
    // (directed only)
    double reciprocity;

    // -- Global networks statistics -------------------------------------------
    /// Global network statistics are calculated during the Kronecker tensor product
    /// process and stored in **first** vertex: access via _g[0].state.__name__
    /// This avoids having to generate the full graph in order to analyse its
    /// properties

    // Network clustering coefficient
    double clustering_global;

    // Network degree sequence
    std::vector<std::pair<size_t, size_t>> degree_sequence;

    // Network degree variance
    double degree_variance;

    // Network diameter
    double diameter;

    // Kronecker generation Hamiltonian (error term)
    double error;

    // Size of largest Kronecker factor
    int largest_comp;

    // Network mean degree
    double mean_degree;

    // Number of Kronecker factors
    size_t num_factors;

    // Number of Pareto points found during grid search
    size_t num_Paretos;

    // Number of vertices in the network
    size_t num_vertices;
};

/// The traits of a vertex are just the traits of a graph entity
using VertexTraits = GraphEntityTraits<VertexState>;

/// A vertex is a graph entity with vertex traits
using Vertex = GraphEntity<VertexTraits>;

/// The vertex container type
using VertexContainer = boost::vecS;

// Edge
struct EdgeState
{
    size_t weight = 1.0;
};

/// The traits of an edge are just the traits of a graph entity
using EdgeTraits = GraphEntityTraits<EdgeState>;

/// An edge is a graph entity with edge traits
using Edge = GraphEntity<EdgeTraits>;

/// The edge container type
using EdgeContainer = boost::vecS;

// Undirected Graph
using NWType = boost::adjacency_list<EdgeContainer,
                                     VertexContainer,
                                     boost::undirectedS,
                                     Vertex,
                                     boost::property<
                                          boost::edge_weight_t,
                                          double,
                                          EdgeState>
                                     >;
}

#endif // UTOPIA_MODELS_KRONGEN_GRAPHDEF
