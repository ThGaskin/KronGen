#define BOOST_TEST_MODULE Kronecker graph utils test

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
    Infrastructure() : BaseInfrastructure<>("test_utils.yml") {};
};

struct VertexState {
    double clustering_global = -1;
    double diameter = -1;
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

BOOST_FIXTURE_TEST_CASE(test_Kronecker_properties, Test_Graph)
{
    using vertices_size_type = typename boost::graph_traits<G_vec_u>::vertices_size_type;

    test_config_callable (

      [&](auto test_cfg){

          G_vec_u g0{};
          auto v0 = boost::add_vertex(g0);
          boost::add_edge(v0, v0, g0);

          std::vector<std::size_t>N;
          std::vector<double>k;
          std::vector<double>var;
          std::vector<double>c;
          std::vector<double>diam;

          BOOST_TEST_CHECKPOINT("Generating Kronecker graph ... ");
          for (const auto& factor_map : test_cfg["Kronecker"]) {
              auto g = create_graph<G_vec_u>(factor_map.second, *rng);

              const auto model = get_as<std::string>("model", factor_map.second);
              const auto n_vertices = boost::num_vertices(g);
              const auto deg_stats = degree_statistics(g);
              const auto clustering = global_clustering_coeff(g);
              const auto starting_point = fourSweep<vertices_size_type>(g);
              const double d = iFUB(starting_point.first, starting_point.second, 0, g);

              N.push_back(n_vertices);
              k.push_back(static_cast<double>(deg_stats.first));
              var.push_back(static_cast<double>(deg_stats.second));
              c.push_back(clustering);
              diam.push_back(d);

              for (const auto w : Utopia::range<IterateOver::vertices>(g)) {
                  add_edge(w, w, g);
              }
              g0 = Kronecker_product(g0, g);

              //Check properties of particular graphs: zero_c, chain
              if (model == "chain") {
                  BOOST_TEST(d == n_vertices - 1);
              }
              else if (model == "zero_c"){
                  BOOST_TEST(clustering == 0);
              }
          }
          for (const auto w : Utopia::range<IterateOver::vertices>(g0)) {
              remove_edge(w, w, g0);
          }

          BOOST_TEST_CHECKPOINT("Kronecker graph generated");

          const auto deg_stats = degree_statistics(g0);

          // Check number of vertices
          const std::size_t n_vertices = num_vertices(g0);
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
          const auto starting_point = fourSweep<vertices_size_type>(g0);
          const auto d = iFUB(starting_point.first, starting_point.second, 0, g0);
          BOOST_TEST_CHECKPOINT("Testing diameter");
          BOOST_TEST(d == std::max(diam[0], diam[1]));

      },
      cfg
    );
}
