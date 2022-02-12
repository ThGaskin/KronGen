#define BOOST_TEST_MODULE Kronecker graph properties test

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include <utopia/core/testtools.hh>
#include <utopia/core/types.hh>
#include <utopia/core/graph.hh>

#include "../aux_graphs.hh"
#include "../KronGen.hh"
#include "../graph_creation.hh"
#include "test_utils.hh"
#include "../utils.hh"

#include "../../NetworkAnalyser/graph_metrics.hh"

using namespace Utopia::TestTools;
using namespace Utopia::Models::KronGen;
using namespace Utopia::Models::KronGen::Utils;
using namespace Utopia::Models::NetworkAnalyser;
using namespace Utopia;

// -- Types --------------------------------------------------------------------

struct Infrastructure : public BaseInfrastructure<> {
    Infrastructure() : BaseInfrastructure<>("test_Kronecker_properties.yml") {};
};

/// The test graph types
struct Test_Graph : Infrastructure {

  // undirected
  using Graph = Utopia::Models::KronGen::NWType;

};

// -- Tests --------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(test_Kronecker_properties, Test_Graph)
{

    test_config_callable (

      [&](auto test_cfg){

          const auto G = GraphCreation::create_graph<Graph>(test_cfg, *rng, log);

          TestUtils::assert_no_parallel_self_edges(G);

          BOOST_TEST_CHECKPOINT("Kronecker graph generated");

          // Check clustering coefficient
          const auto c0 = global_clustering_coeff(G);
          BOOST_TEST_CHECKPOINT("Testing clustering coefficient");
          BOOST_TEST(c0 == G[0].state.clustering_global,
                     boost::test_tools::tolerance(1.e-12));

          // Check degree sequence
          const auto deg_seq = degree_sequence(G);
          BOOST_TEST_CHECKPOINT("Testing degree sequence");
          BOOST_TEST(deg_seq == G[0].state.degree_sequence);

          const auto deg_stats = degree_statistics(G);

          // Check degree distribution variance
          const double variance = deg_stats.second;
          BOOST_TEST_CHECKPOINT("Testing degree distribution variance");
          BOOST_TEST(variance == G[0].state.degree_variance,
                     boost::test_tools::tolerance(1.e-12));

          // Check diameter
          const auto d = diameter(G);
          BOOST_TEST_CHECKPOINT("Testing diameter");
          BOOST_TEST(d == G[0].state.diameter);

          // Check mean degree
          const double mean_degree = deg_stats.first;
          BOOST_TEST_CHECKPOINT("Testing mean degree");
          BOOST_TEST(mean_degree == G[0].state.mean_degree,
                     boost::test_tools::tolerance(1.e-12));

          // Check number of vertices
          const std::size_t n_vertices = boost::num_vertices(G);
          BOOST_TEST_CHECKPOINT("Testing vertex count");
          BOOST_TEST(n_vertices == G[0].state.num_vertices);

      },
      cfg
    );
}
