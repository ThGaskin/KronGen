#define BOOST_TEST_MODULE Kronecker graph utils test

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

using namespace Utopia::Models::KronGen::GraphCreation;
using namespace Utopia::Models::KronGen::Utils;
using namespace Utopia::Models::NetworkAnalyser;
using namespace Utopia;
using namespace Utopia::TestTools;

/* Missing: diameter estimation
*/

// -- Types --------------------------------------------------------------------
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
  using Graph = boost::adjacency_list<
                      boost::vecS,         // edge container
                      boost::vecS,         // vertex container
                      boost::undirectedS,
                      Vertex,              // vertex struct
                      Edge>;               // edge struct

};

// -- Tests --------------------------------------------------------------------
// Test factorisation of N (number of vertices)
BOOST_AUTO_TEST_CASE (test_N_factors)
{
    std::vector<std::size_t> N = {15, 16, 20, 25, 30, 35, 40, 100};
    std::vector<std::size_t> primes = {0, 1, 2, 3, 4, 5, 7, 10, 11, 14, 17, 19, 22, 167};

    for (const auto& n : N) {
        const auto factors = N_factors(n);
        const auto closest_factors = closest_N_factors(n);
        BOOST_TEST(!factors.empty());
        BOOST_TEST(factors.size() == closest_factors.size());
        for (std::size_t i = 0; i < factors.size(); ++i) {
            const auto f = factors[i];
            BOOST_TEST(f.first*f.second == n);
            BOOST_TEST(f.first != 1);
            BOOST_TEST(f.first != n);
            BOOST_TEST(f.first == closest_factors[i].first);
            BOOST_TEST(f.second == closest_factors[i].second);
        }
    }

    for (const auto& p : primes) {
        const auto factors = N_factors(p);
        BOOST_TEST(factors.empty());
        const auto closest_factors = closest_N_factors(p);
        BOOST_TEST(!closest_factors.empty());
        for (const auto& f : closest_factors) {
            if (p < 9) {
                BOOST_TEST(!(9 % f.first));
                BOOST_TEST(!(9 % f.second));
            }
            else {
                BOOST_TEST((!((p-1) % f.first) or !((p+1) % f.first)));
                BOOST_TEST((!((p-1) % f.second) or !((p+1) % f.second)));
            }
        }
    }
}

// Test clustering coefficient and diameter are analysed and written into graph
// vertices
BOOST_FIXTURE_TEST_CASE(test_Kronecker_analysis, Test_Graph)
{

    test_config_callable (

      [&](auto test_cfg){

          const Graph g = create_graph<Graph>(test_cfg, *rng, logger, true);

          BOOST_TEST_CHECKPOINT("Kronecker graph generated");

          // Check clustering
          const auto c0 = global_clustering_coeff(g);
          BOOST_TEST_CHECKPOINT("Testing clustering coefficient");
          BOOST_TEST(g[0].state.clustering_global == c0,
                     boost::test_tools::tolerance(1.e-12)
          );

          // Check diameter
          const auto d = diameter(g);
          BOOST_TEST_CHECKPOINT("Testing diameter");
          BOOST_TEST(g[0].state.diameter == d);
      },
      cfg
    );
}

// Test factorisation of m (mean degree)
BOOST_AUTO_TEST_CASE (test_mean_deg_factors)
{
    std::vector<std::size_t> degrees = {0, 1, 2, 3, 4, 5, 7, 9, 10, 12, 14, 16, 19, 20, 21, 50};

    for (const auto& m : degrees) {
        const auto factors = mean_deg_factors(m);
        const auto closest_factors = closest_mean_deg_factors(m);
        if (m < 3) {
            BOOST_TEST(factors.empty());
            BOOST_TEST(!closest_factors.empty());
            for (const auto& f : closest_factors) {
                BOOST_TEST(f.first == 2);
                BOOST_TEST(f.second == 2);
            }
        }
        if (factors.empty()) {
            BOOST_TEST(!(closest_factors.empty()));
        }
        else {
            BOOST_TEST(closest_factors.size() == factors.size());
        }
        for (const auto& f : closest_factors) {
            const auto mm = Kronecker_mean_degree(f.first, f.second);
            BOOST_TEST(Kronecker_mean_degree_inv(mm, f.first) == f.second);
            BOOST_TEST(Kronecker_mean_degree_inv(mm, f.second) == f.first);
        }
    }
}

// Test clustering coefficient of regular graph
BOOST_FIXTURE_TEST_CASE (test_regular_graph_properties, Test_Graph)
{
    const size_t N = 200;
    std::vector<std::size_t> degrees = {4, 6, 8, 10, 20, 40, 100, 180};
    for (const auto& m : degrees) {
        const Graph g = Utopia::Graph::create_regular_graph<Graph>(N, m, false);
        const auto c = global_clustering_coeff(g);
        BOOST_TEST(c == regular_graph_clustering(N, m),
                   boost::test_tools::tolerance(0.05));
    }
}
