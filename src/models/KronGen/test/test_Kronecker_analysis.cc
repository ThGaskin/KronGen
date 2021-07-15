#define BOOST_TEST_MODULE Kronecker graph analysis test

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/clustering_coefficient.hpp>
#include <spdlog/spdlog.h>

#include <utopia/core/testtools.hh>
#include <utopia/core/types.hh>
#include <utopia/core/graph.hh>

#include "../graph_creation.hh"
#include "../../NetworkAnalyser/graph_metrics.hh"
#include "../utils.hh"

using namespace Utopia::TestTools;
using namespace Utopia::Models::KronGen::GraphCreation;
using namespace Utopia::Models::KronGen::Utils;
using namespace Utopia::Models::NetworkAnalyser;
using namespace Utopia;

// -- Types -------------------------------------------------------------------
auto logger = Utopia::init_logger("root.KronGen", spdlog::level::info);

struct Infrastructure : public BaseInfrastructure<> {
    Infrastructure() : BaseInfrastructure<>("test_Kronecker_analysis.yml") {};
};

struct VertexState {
    double clustering_global = -1;
    double diameter = -1;
    double var = -1;
    double mean_deg = -1;
};
struct Edge {};

/// The test graph types
struct Test_Graph : Infrastructure {

  using VertexTraits = Utopia::GraphEntityTraits<VertexState>;
  using Vertex = Utopia::GraphEntity<VertexTraits>;

  // undirected
  using G_vec_u = boost::adjacency_list<
                      boost::vecS,         // edge container
                      boost::vecS,         // vertex container
                      boost::undirectedS,
                      Vertex,              // vertex struct
                      Edge>;               // edge struct

};
