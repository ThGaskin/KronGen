#ifndef UTOPIA_MODELS_KRONGEN_TYPEDEFS
#define UTOPIA_MODELS_KRONGEN_TYPEDEFS

#include <any>
#include <map>
#include <string>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include <utopia/core/model.hh>
#include <utopia/core/types.hh>
#include <utopia/core/graph/entity.hh>

// .. Type definitions used in the model .......................................

namespace Utopia::Models::KronGen::TypeDefinitions{

// A factor is a list of ints
using factor = typename std::vector<size_t>;

// Factors are lists of factors
using factors = typename std::vector<factor>;

// Vectors of pairs of integers
using vector_pt = typename std::vector<std::pair<size_t, size_t>>;

// Entry type of the map of analysis targets: an entry is a map consisting of
// strings as keys ("target", "calculate", etc.) and a bool and value pair.
// The bool indicates whether or not a parameter is a target/analysis entry,
// and the value indicates the target value/calculated analysis value.
// Target values are const.
using entry_type = typename std::map<std::string, std::pair<bool, std::any>>;

// A map consists of key and entry_type pairs
using map_type = typename std::map<std::string, entry_type>;

// All graph types considered in the KronGen model
enum GraphType {
    Chain,
    Complete,
    ErdosRenyi,
    KlemmEguiluz,
    Regular,
    BarabasiAlbert,
    SmallWorld
};

// Simple struct containing the information required to generate a graph
struct GraphDesc {

    size_t num_vertices;

    size_t mean_degree;

    GraphType type;

    double mu;
};

// A Pareto point is a list of graph descriptors
using ParetoPoint = typename std::vector<GraphDesc>;

// A Pareto set if a list of Pareto points
using ParetoSet = typename std::vector<ParetoPoint>;

// Convenient printing function for graph types
std::string Graph_Type[] = {
    "Chain",
    "Complete",
    "Erdos-Renyi",
    "Klemm-Eguiluz",
    "Regular",
    "Barabasi-Albert",
    "Small-world"
};

// Converts a string to a GraphType
GraphType to_graphtype(std::string s) {
    if (s == "BarabasiAlbert"){
        return BarabasiAlbert;
    }
    else if (s == "chain") {
        return Chain;
    }
    else if (s == "complete") {
        return Complete;
    }
    else if (s == "ErdosRenyi") {
        return ErdosRenyi;
    }
    else if (s == "KlemmEguiluz") {
        return KlemmEguiluz;
    }
    else if (s == "regular" or s == "Regular") {
        return Regular;
    }
    else if (s == "WattsStrogatz"){
        return SmallWorld;
    }
    else {
        return ErdosRenyi;
    }
}

// The graph vertex type
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

#endif // UTOPIA_MODELS_KRONGEN_TYPEDEFS
