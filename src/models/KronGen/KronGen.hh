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
#include "graph_definition.hh"

// The NetworkAnalyser
#include "../NetworkAnalyser/NetworkAnalyser.hh"

namespace Utopia::Models::KronGen {

using namespace Utopia::Models::KronGen::GraphDefinition;

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

using vector = typename std::vector<double>;

/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<>;

// Network analyser
using NWAnalyser = NetworkAnalyser::NetworkAnalyser<NWType>;

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
        typename boost::graph_traits<NWType>::vertex_descriptor;

    /// Data type for an edge descriptor
    using EdgeDesc = typename boost::graph_traits<NWType>::edge_descriptor;

    /// Data type for a rule function operating on vertices returning void
    using VertexVoidRule = typename std::function<void(VertexDesc, NWType&)>;

    /// Data type for a rule function operating on vertices returning a state
    using VertexStateRule =
        typename std::function<VertexState(VertexDesc, NWType&)>;

    /// Data type for a rule function operating on edges returning void
    using EdgeVoidRule = typename std::function<void(EdgeDesc, NWType&)>;

    /// Data type for a rule function operating on edges returning a state
    using EdgeStateRule =
        typename std::function<EdgeState(EdgeDesc, NWType&)>;

private:

  NWType _g;

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

    NWType initialize_graph() {

      this->_log->info("Creating the graph ...");

      NWType g = GraphCreation::create_graph<NWType>(
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
