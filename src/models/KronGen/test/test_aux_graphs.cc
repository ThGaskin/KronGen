#define BOOST_TEST_MODULE Kronecker aux graph test

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/clustering_coefficient.hpp>

#include <utopia/core/testtools.hh>
#include <utopia/core/types.hh>
#include <utopia/core/graph.hh>

#include "../graph_creation.hh"
#include "../graph_properties.hh"
#include "../KronGen.hh"
#include "../../NetworkAnalyser/graph_metrics.hh"
#include "test_utils.hh"
#include "../utils.hh"

using namespace boost;
using namespace Utopia::TestTools;
using namespace Utopia::Models::KronGen::AuxGraphs;
using namespace Utopia::Models::KronGen::TestUtils;
using namespace Utopia::Models::KronGen::Utils;
using namespace Utopia::Models::KronGen::GraphProperties;
using namespace Utopia::Models::NetworkAnalyser;
using namespace Utopia;

auto logger = Utopia::init_logger("root.KronGen", spdlog::level::info);

struct Infrastructure : public BaseInfrastructure<> {
    Infrastructure() : BaseInfrastructure<>("test_aux_graphs.yml") {};
};

/// The test graph types
struct Test_Graph : Infrastructure {

  // undirected
  using Graph = Utopia::Models::KronGen::NWType;

};

// Creation of the zero clustering graph, a bipartite k-regular graph with N vertices
BOOST_FIXTURE_TEST_CASE(zero_c_graph, Test_Graph)
{
    const std::vector<double> mean_degree = {2, 3, 4, 5, 6, 8, 10, 11, 20};
    for (const auto& k : mean_degree) {
        for (std::size_t N = 2*k; N < 4*k; N+=2){
            const auto G = create_zero_c_graph<Graph>(N, k);
            const auto c = global_clustering_coeff(G);

            BOOST_TEST(num_vertices(G) == N);
            BOOST_TEST(2*num_edges(G) == k*N);
            BOOST_TEST(c == 0);
            assert_no_parallel_self_edges(G);
        }
    }
}

// Creation of a chain graph of length N
BOOST_FIXTURE_TEST_CASE(chain_graph, Test_Graph)
{
    const std::vector<std::size_t> n_vertices = {2, 3, 4, 5, 6, 8, 10, 11, 20};
    for (const auto& N : n_vertices) {
        const auto g = create_chain_graph<Graph>(N);

        BOOST_TEST(num_vertices(g) == N);
        BOOST_TEST(2.0*num_edges(g)/N == mean_degree_chain(N));
        BOOST_TEST(diameter(g) == N - 1);
        assert_no_parallel_self_edges(g);
    }
}

// Creation of a star graph
BOOST_FIXTURE_TEST_CASE(star_graph, Test_Graph)
{
    const std::vector<std::size_t> n_vertices = {20, 30, 40};
    const std::vector<double> mean_degree = {2, 3, 4, 5, 6, 7, 8};
    const std::vector<std::size_t> diameters = {2, 3, 4, 5, 6, 7, 8};

    for (const auto& N : n_vertices) {
        for (const auto& k : mean_degree) {
            for (const auto& d : diameters) {
                const auto g = create_star_graph<Graph>(N, k, d, *rng);

                BOOST_TEST(num_vertices(g) == N);
                BOOST_TEST(2*num_edges(g) == N*k);
                BOOST_TEST(d == diameter(g));
                assert_no_parallel_self_edges(g);
            }
        }
    }
}

// Test config extraction
BOOST_FIXTURE_TEST_CASE(graph_extraction, Test_Graph)
{
  test_config_callable (

    [&](auto test_cfg){

      const auto G = Utopia::Models::KronGen::AuxGraphs::create_graph<Graph>(test_cfg, *rng);
      const size_t N = get_as<size_t>("num_vertices", test_cfg);
      BOOST_TEST(boost::num_vertices(G) == N);

    },
    cfg
  );
}
