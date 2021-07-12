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
#include "../../NetworkAnalyser/graph_metrics.hh"
#include "test_utils.hh"
#include "../utils.hh"

using namespace boost;
using namespace Utopia::TestTools;
using namespace Utopia::Models::KronGen::AuxGraphs;
using namespace Utopia::Models::KronGen::TestUtils;
using namespace Utopia::Models::KronGen::Utils;
using namespace Utopia::Models::NetworkAnalyser;
using namespace Utopia;

struct Infrastructure : public BaseInfrastructure<> {
    Infrastructure() : BaseInfrastructure<>() {};
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
  using Graph = boost::adjacency_list<
                      boost::vecS,         // edge container
                      boost::vecS,         // vertex container
                      boost::undirectedS,
                      Vertex,              // vertex struct
                      Edge>;               // edge struct

  using vertices_size_type = typename boost::graph_traits<Graph>::vertices_size_type;

};

// zero clustering graph
BOOST_FIXTURE_TEST_CASE(zero_c_graph, Test_Graph)
{
    const std::vector<double> mean_degree = {2, 3, 4, 5, 6, 8, 10, 11, 20};
    for (const auto& k : mean_degree) {
        for (std::size_t N = 2*k; N < 4*k; N+=2){
            const auto g = create_zero_c_graph<Graph>(N, k);
            const auto c = global_clustering_coeff(g);

            BOOST_TEST(num_vertices(g) == N);
            BOOST_TEST(2*num_edges(g) == k*N);
            BOOST_TEST(c==0);
            assert_no_parallel_self_edges(g);
        }
    }
}

// Chain graph
BOOST_FIXTURE_TEST_CASE(chain_graph, Test_Graph)
{
    const std::vector<std::size_t> n_vertices = {2, 3, 4, 5, 6, 8, 10, 11, 20};
    for (const auto& N : n_vertices) {
        const auto g = create_chain_graph<Graph>(N);
        const auto starting_point = fourSweep<vertices_size_type>(g);
        const auto diam = iFUB(starting_point.first, starting_point.second, 0, g);

        BOOST_TEST(num_vertices(g) == N);
        BOOST_TEST(num_edges(g) == num_vertices(g)-1);
        BOOST_TEST(diam == N - 1);
        assert_no_parallel_self_edges(g);
    }
}

// Star graph
BOOST_FIXTURE_TEST_CASE(star_graph, Test_Graph)
{

    const std::vector<std::size_t> n_vertices = {20, 30, 40};
    const std::vector<double> mean_degree = {2, 3, 4, 5, 6, 7, 8};
    const std::vector<std::size_t> diameter = {2, 3, 4, 5, 6, 7, 8};

    for (const auto& N : n_vertices) {
        for (const auto& k : mean_degree) {
            for (const auto& d : diameter) {
                const auto g = create_star_graph<Graph>(N, k, d, *rng);
                const auto starting_point = fourSweep<vertices_size_type>(g);
                const auto diam = iFUB(starting_point.first, starting_point.second, 0, g);

                BOOST_TEST(num_vertices(g) == N);
                BOOST_TEST(2*num_edges(g) == N*k);
                BOOST_TEST(diam == d);
                assert_no_parallel_self_edges(g);
            }
        }
    }
}
