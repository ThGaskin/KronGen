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

#include "../KronGen.hh"
#include "../graph_creation.hh"
#include "../graph_properties.hh"
#include "../../NetworkAnalyser/graph_metrics.hh"
#include "../utils.hh"

using namespace Utopia;
using namespace Utopia::Models::KronGen::GraphCreation;
using namespace Utopia::Models::KronGen::GraphProperties;
using namespace Utopia::Models::KronGen::Utils;
using namespace Utopia::Models::NetworkAnalyser;
using namespace Utopia::TestTools;


// -- Types --------------------------------------------------------------------
auto logger = Utopia::init_logger("root.KronGen", spdlog::level::info);

struct Infrastructure : public BaseInfrastructure<> {
    Infrastructure() : BaseInfrastructure<>("test_Kronecker_analysis.yml") {};
};

/// The test graph types
struct Test_Graph : Infrastructure {

  // undirected
  using Graph = Utopia::Models::KronGen::GraphType;

};

// -- Tests --------------------------------------------------------------------
// Test getting the grid in N
BOOST_AUTO_TEST_CASE (test_N_grid)
{
    const std::vector<size_t> targets = {1, 2, 169, 500, 1000};
    const std::vector<double> errors = {0., 0.1, 0.2};
    const size_t min_dim = 2;
    const size_t max_dim = 8;
    for (const auto& N_target : targets){
        for (const auto& err : errors) {

            const auto N_grid = get_N_grid (N_target, err, min_dim, max_dim);

            for (const auto& n : N_grid) {
                size_t N_res = 1;
                for (const auto& fac : n) {
                    N_res *= fac;
                }
                BOOST_TEST(1.0*N_res >= (N_target*(1.-err)));
                BOOST_TEST(1.0*N_res <= (N_target*(1.+err)));
                BOOST_TEST(n.size() >= min_dim);
                BOOST_TEST(n.size() <= max_dim);
            }
        }
    }
}


// Test getting the grid in k
BOOST_AUTO_TEST_CASE (test_k_grid)
{
    const std::vector<size_t> targets = {1, 2, 51, 169, 500};
    const std::vector<double> errors = {0., 0.1, 0.2};
    const size_t min_dim = 2;
    const size_t max_dim = 8;

    for (const auto& k_target : targets){
        for (const auto& err : errors) {

            const auto k_grid = get_k_grid (k_target, err, min_dim, max_dim);

            for (const auto& k : k_grid) {
                size_t k_res = 1;
                for (const auto& fac : k) {
                    k_res *= (fac+1);
                }
                k_res -=1;
                BOOST_TEST(1.0*k_res >= (k_target*(1.-err)));
                BOOST_TEST(1.0*k_res <= (k_target*(1.+err)));
                BOOST_TEST(k.size() >= min_dim);
                BOOST_TEST(k.size() <= max_dim);
                if (k.size() == 2) {
                    BOOST_TEST(Kronecker_mean_degree_inv(k_res, k[1]) == k[0]);
                }
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

// Test properties of regular graph
BOOST_FIXTURE_TEST_CASE (test_regular_graph_properties, Test_Graph)
{
    const size_t N = 200;
    std::vector<std::size_t> degrees = {4, 6, 8, 10, 20, 40, 100, 180};
    for (const auto& m : degrees) {
        const Graph g = Utopia::Graph::create_regular_graph<Graph>(N, m, false);
        const auto c = global_clustering_coeff(g);
        BOOST_TEST(c == clustering_regular(N, m),
                   boost::test_tools::tolerance(0.05));
    }
}

// // Test Erdos-Renyi properties
// BOOST_FIXTURE_TEST_CASE(create_ER_graph, Test_Graph)
// {
//     const double N = 300;
//     const std::vector<double> mean_deg = {10, 20, 80};
//
//     for (const auto& k : mean_deg){
//         const auto g = Utopia::Graph::create_ErdosRenyi_graph<Graph>(N, k, false, false, *rng);
//         const auto var = degree_statistics(g).second;
//         const auto c = global_clustering_coeff(g);
//         BOOST_TEST(var == degree_variance_ER(N, k), boost::test_tools::tolerance(0.1));
//         BOOST_TEST(c == clustering_ER(N, k), boost::test_tools::tolerance(0.12));
//     }
// }

// Test objective functions
