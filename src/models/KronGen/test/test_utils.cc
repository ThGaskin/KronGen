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
using namespace Utopia;

// -- Types -------------------------------------------------------------------

struct Infrastructure : public BaseInfrastructure<> {
    Infrastructure() : BaseInfrastructure<>("test_utils.yml") {};
};

struct Vertex {};
struct Edge {};

/// The test graph types
struct Test_Graph : Infrastructure {

  // undirected
  using G_vec_u = boost::adjacency_list<
                      boost::vecS,         // edge container
                      boost::vecS,         // vertex container
                      boost::undirectedS,
                      Vertex,              // vertex struct
                      Edge>;               // edge struct
};

BOOST_FIXTURE_TEST_CASE(test_clustering_coeff, Test_Graph)
{
    test_config_callable (

      [&](auto test_cfg){

          G_vec_u g0{};
          auto v0 = boost::add_vertex(g0);
          boost::add_edge(v0, v0, g0);

          std::vector<double>c;
          std::vector<double>k;
          std::vector<double>v;

          for (const auto& factor_map : test_cfg["Kronecker"]) {
              auto g = create_graph<G_vec_u>(factor_map.second, *rng);
              const auto deg_stats = Utopia::Models::NetworkAnalyser::degree_statistics(g);
              c.push_back(Utopia::Models::NetworkAnalyser::global_clustering_coeff(g));
              k.push_back(static_cast<double>(deg_stats.first));
              v.push_back(static_cast<double>(deg_stats.second));
              for (const auto w : Utopia::range<IterateOver::vertices>(g)) {
                  add_edge(w, w, g);
              }
              g0 = Kronecker_product(g0, g);
          }
          for (const auto w : Utopia::range<IterateOver::vertices>(g0)) {
              remove_edge(w, w, g0);
          }

          const auto c_r = Kronecker_clustering(c[0], c[1], k[0], k[1], v[0], v[1]);
          const auto c0 = Utopia::Models::NetworkAnalyser::global_clustering_coeff(g0);

          BOOST_TEST(c_r == c0, boost::test_tools::tolerance(1.e-12));

      },
      cfg
    );
}
