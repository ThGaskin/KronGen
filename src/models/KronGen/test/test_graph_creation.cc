#define BOOST_TEST_MODULE Kronecker graph creation test

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/clustering_coefficient.hpp>

#include <utopia/core/testtools.hh>
#include <utopia/core/types.hh>
#include <utopia/core/graph.hh>

#include "../graph_creation.hh"

using namespace Utopia::Graph;
using namespace Utopia::TestTools;
using namespace Utopia::Models::KronGen::GraphCreation;


// -- Types -------------------------------------------------------------------

struct Infrastructure : public BaseInfrastructure<> {
    Infrastructure() : BaseInfrastructure<>("test_graph_creation.yml") {};
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

// -- Helper functions --------------------------------------------------------
// Test against parallel or self edges
template<typename Graph>
void assert_no_parallel_self_edges(Graph& G) {
    size_t num_parallel = 0;
    for (auto [v, v_end] = boost::vertices(G); v!=v_end; ++v) {
        for (auto [e, e_end] = boost::out_edges(*v, G); e!=e_end; ++e) {
            int counter = 0;
            for (auto [g, g_end] = boost::out_edges(*v, G); g!=g_end; ++g) {
                if (target(*g, G) == target(*e, G)) {
                    counter += 1;
                }
            }
            if (counter > 1) {
                num_parallel += 1;
            }
            // Check against self-edges
            BOOST_TEST(target(*e, G) != *v);
        }
    }
    BOOST_TEST(num_parallel == 0);
}

// -- Tests -------------------------------------------------------------------
// Undirected graphs
BOOST_FIXTURE_TEST_CASE(create_Kron_graph, Test_Graph)
{
    test_config_callable (

      [&](auto test_cfg){

        const auto g0 = create_Kronecker_graph<G_vec_u>(test_cfg, *rng);
        std::size_t num_vertices = 1;
        std::size_t mean_degree = 1;
        bool is_complete = true;

        for (const auto& factor_map : test_cfg["Kronecker"]) {

            const auto factor_cfg = factor_map.second;

            num_vertices *= get_as<std::size_t>("num_vertices", factor_cfg);
            if (get_as<std::string>("model", factor_cfg) == "complete") {
                mean_degree *= get_as<std::size_t>("num_vertices", factor_cfg);
            }
            else {
                mean_degree *= (get_as<std::size_t>("mean_degree", factor_cfg)+1);
                is_complete = false;
            }
        }
        mean_degree -= 1;
        const auto num_edges = std::round(static_cast<double>(num_vertices)/2*mean_degree);
        BOOST_TEST(boost::num_vertices(g0) == num_vertices);
        BOOST_TEST(boost::num_edges(g0) == num_edges);

        assert_no_parallel_self_edges(g0);

        if (is_complete) {
            for (auto [v, v_end] = boost::vertices(g0); v!=v_end; ++v) {
                double c = boost::clustering_coefficient(g0, *v);
                BOOST_TEST(c == 1);
            }
        }

      },
      cfg
    );
}

BOOST_FIXTURE_TEST_CASE(create_zero_clustering_graph, Test_Graph)
{
    const std::vector<std::size_t> n_vertices{4, 6, 8, 10, 12, 14, 16};

    for (const auto& N : n_vertices) {
        for (std::size_t k = 2; k < (N/2+1); ++k ){
            const auto g = create_zero_c_graph<G_vec_u>(N, k);
            double deg_sum = 0;
            BOOST_TEST(boost::num_vertices(g) == N);
            for (auto [v, v_end] = boost::vertices(g); v!=v_end; ++v) {
                double c = boost::clustering_coefficient(g, *v);
                auto deg = degree(*v, g);
                BOOST_TEST(c == 0);
                deg_sum += deg;
                BOOST_TEST(deg == k);
            }
            BOOST_TEST(deg_sum/N == k);
            BOOST_TEST(Utopia::Models::NetworkAnalyser::global_clustering_coeff(g) == 0);

            assert_no_parallel_self_edges(g);
        }
    }
}
