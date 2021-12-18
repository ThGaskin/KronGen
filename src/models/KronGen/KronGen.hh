#ifndef UTOPIA_MODELS_KRONGEN_HH
#define UTOPIA_MODELS_KRONGEN_HH

// standard library includes
#include <random>

// third-party library includes
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>

// Kronecker Graph generation file
#include "graph_creation.hh"

// The NetworkAnalyser
#include "../NetworkAnalyser/NetworkAnalyser.hh"

namespace Utopia::Models::KronGen {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

using vector = typename std::vector<double>;

// -- Vertex ------------------------------------------------------------------

struct VertexState
{
    // Betweenness centrality
    double betweenness;

    // Closeness centrality
    double closeness;

    // The vertex' clustering coefficient, for not having to compute repeatedly
    double clustering_coeff;

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

    // global statistics: can be calculated during Kronecker tensor product
    // process and stored in **first** vertex: access via _g[0].state.__name__
    double clustering_global = -1;
    double diameter = -1;
    double error; // optimisation error
    int largest_comp = -1;
    double mean_deg = -1;
    size_t num_factors = 1; // number of Kronecker factors
    size_t num_Paretos = 1;
    int num_vertices = -1;
    double var = -1;

};

/// The traits of a vertex are just the traits of a graph entity
using VertexTraits = GraphEntityTraits<VertexState>;

/// A vertex is a graph entity with vertex traits
using Vertex = GraphEntity<VertexTraits>;

/// The vertex container type
using VertexContainer = boost::vecS;

// -- Edge --------------------------------------------------------------------

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

// -- Graph -------------------------------------------------------------------
// undirected graph

using GraphType = boost::adjacency_list<EdgeContainer,
                                        VertexContainer,
                                        boost::undirectedS,
                                        Vertex,
                                        boost::property<
                                          boost::edge_weight_t,
                                          double,
                                          EdgeState>>;

/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<>;

// Network analyser
using NWAnalyser = NetworkAnalyser::NetworkAnalyser<GraphType>;

// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// The KronGen Model
class KronGen: public Model<KronGen, ModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<KronGen, ModelTypes>;

    /// Data type of the group to write model data to, holding datasets
    using DataGroup = typename Base::DataGroup;

    /// Data type for a dataset
    using DataSet = typename Base::DataSet;

    // .. Graph-related types and rule types ..................................

    /// Data type for a vertex descriptor
    using VertexDesc =
        typename boost::graph_traits<GraphType>::vertex_descriptor;

    /// Data type for an edge descriptor
    using EdgeDesc = typename boost::graph_traits<GraphType>::edge_descriptor;

    /// Data type for a rule function operating on vertices returning void
    using VertexVoidRule = typename std::function<void(VertexDesc, GraphType&)>;

    /// Data type for a rule function operating on vertices returning a state
    using VertexStateRule =
        typename std::function<VertexState(VertexDesc, GraphType&)>;

    /// Data type for a rule function operating on edges returning void
    using EdgeVoidRule = typename std::function<void(EdgeDesc, GraphType&)>;

    /// Data type for a rule function operating on edges returning a state
    using EdgeStateRule =
        typename std::function<EdgeState(EdgeDesc, GraphType&)>;

private:

  GraphType _g;

  NWAnalyser _nwanalyser;


public:
    // -- Model Setup ---------------------------------------------------------

    /// Construct the KronGen model
    template<class ParentModel>
    KronGen (
        const std::string& name,
        ParentModel& parent_model,
        const DataIO::Config& custom_cfg = {}
    )
    :
        Base(name, parent_model, custom_cfg),

        _g(initialize_graph()),

        _nwanalyser("NetworkAnalyser", *this, _g)

    {}


private:
    // .. Setup functions .....................................................

    GraphType initialize_graph() {

      this->_log->info("Creating the graph ...");

      GraphType g = GraphCreation::create_graph<GraphType>(
            this->_cfg,
            *this->_rng,
            this->_log);

      this->_log->info("Graph creation complete.");

      return g;
    }

    // .. Helper functions ....................................................

public:
    // -- Public Interface ----------------------------------------------------

    void prolog () {
        _nwanalyser.prolog();
    }
    void perform_step () {}
    void monitor () {}
    void write_data () {}
    void epilog () {
        _nwanalyser.epilog();
    }

    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

};

} // namespace Utopia::Models::KronGen

#endif // UTOPIA_MODELS_KRONGEN_HH
