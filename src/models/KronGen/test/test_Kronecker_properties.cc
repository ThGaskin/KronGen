#define BOOST_TEST_MODULE Kronecker graph properties test

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include <utopia/core/testtools.hh>
#include <utopia/core/types.hh>
#include <utopia/core/graph.hh>

#include "../aux_graphs.hh"
#include "../graph_creation.hh"
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
  using Graph = boost::adjacency_list<
                      boost::vecS,         // edge container
                      boost::vecS,         // vertex container
                      boost::undirectedS,
                      Vertex,              // vertex struct
                      Edge>;               // edge struct

};

// -- Tests --------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(test_Kronecker_properties, Test_Graph)
{
  
    test_config_callable (

      [&](auto test_cfg){

          Graph g0{1};
          add_self_edges(g0);

          std::vector<std::size_t>N;
          std::vector<double>k;
          std::vector<double>var;
          std::vector<double>c;
          std::vector<double>diam;

          BOOST_TEST_CHECKPOINT("Generating Kronecker graph ... ");
          for (const auto& factor_map : test_cfg["Kronecker"]) {
              auto g = GraphCreation::create_graph<Graph>(factor_map.second, *rng, log);

              const auto model = get_as<std::string>("model", factor_map.second);
              const auto n_vertices = boost::num_vertices(g);
              const auto deg_stats = degree_statistics(g);
              const auto clustering = global_clustering_coeff(g);
              const double d = diameter(g);

              N.push_back(n_vertices);
              k.push_back(static_cast<double>(deg_stats.first));
              var.push_back(static_cast<double>(deg_stats.second));
              c.push_back(clustering);
              diam.push_back(d);

              add_self_edges(g);
              g0 = Kronecker_product(g0, g);

          }
          remove_self_edges(g0);

          BOOST_TEST_CHECKPOINT("Kronecker graph generated");

          const auto deg_stats = degree_statistics(g0);

          // Check number of vertices
          const std::size_t n_vertices = boost::num_vertices(g0);
          BOOST_TEST_CHECKPOINT("Testing vertex count");
          BOOST_TEST(n_vertices == Kronecker_num_vertices(N[0], N[1]));

          // Check mean degree
          const double mean_degree = deg_stats.first;
          BOOST_TEST_CHECKPOINT("Testing mean degree");
          BOOST_TEST(mean_degree == Kronecker_mean_degree(k[0], k[1]),
                     boost::test_tools::tolerance(1.e-12));

          // Check degree distribution variance
          const double variance = deg_stats.second;
          BOOST_TEST_CHECKPOINT("Testing degree distribution variance");
          BOOST_TEST(variance == Kronecker_degree_variance(k[0], k[1], var[0], var[1]),
                     boost::test_tools::tolerance(1.e-12));

          // Check clustering coefficient
          const auto c_t = Kronecker_clustering(c[0], c[1], k[0], k[1], var[0], var[1]);
          const auto c0 = global_clustering_coeff(g0);
          BOOST_TEST_CHECKPOINT("Testing clustering coefficient");
          BOOST_TEST(c_t == c0, boost::test_tools::tolerance(1.e-12));

          // Check inversion function of clustering
          const auto g = get_mean_deg_c(c[0], c_t, k[0], var[0]);
          const auto c_G = c_t >= c[0] ? 1 : 0;
          BOOST_TEST_CHECKPOINT("Testing clustering inversion function");
          BOOST_TEST(Kronecker_clustering(c[0], c_G, k[0], g, var[0], 0) == c_t,
                     boost::test_tools::tolerance(1.e-12));

          // Check diameter
          const auto d = diameter(g0);
          BOOST_TEST_CHECKPOINT("Testing diameter");
          BOOST_TEST(d == std::max(diam[0], diam[1]));

      },
      cfg
    );
}
