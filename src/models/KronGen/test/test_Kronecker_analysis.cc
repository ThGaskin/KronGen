#define BOOST_TEST_MODULE Kronecker graph analysis test

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/clustering_coefficient.hpp>

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

BOOST_FIXTURE_TEST_CASE(test_Kronecker_analysis, Test_Graph)
{
    using vertices_size_type = typename boost::graph_traits<G_vec_u>::vertices_size_type;

    test_config_callable (

      [&](auto test_cfg){

          const G_vec_u g = create_graph<G_vec_u>(test_cfg, *rng, true);

          BOOST_TEST_CHECKPOINT("Kronecker graph generated");

          // Check clustering
          const auto c0 = global_clustering_coeff(g);
          BOOST_TEST_CHECKPOINT("Testing clustering coefficient");
          BOOST_TEST(g[0].state.clustering_global == c0,
                     boost::test_tools::tolerance(1.e-12)
          );

          // Check diameter
          const auto starting_point = fourSweep<vertices_size_type>(g);
          const auto d = iFUB(starting_point.first, starting_point.second, 0, g);
          BOOST_TEST_CHECKPOINT("Testing diameter");
          BOOST_TEST(g[0].state.diameter == d);
      },
      cfg
    );
}
