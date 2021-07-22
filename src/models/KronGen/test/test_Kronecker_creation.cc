#define BOOST_TEST_MODULE Kronecker graph creation test

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/clustering_coefficient.hpp>
#include <spdlog/spdlog.h>


#include <utopia/core/graph.hh>
#include "utopia/core/graph/iterator.hh"
#include <utopia/core/testtools.hh>
#include <utopia/core/types.hh>

#include "../graph_creation.hh"
#include "test_utils.hh"

using namespace Utopia::TestTools;
using namespace Utopia::Models::KronGen;


// -- Types -------------------------------------------------------------------
// Initialise the core logger
auto logger = Utopia::init_logger("root.KronGen", spdlog::level::info);

struct Infrastructure : public BaseInfrastructure<> {
    Infrastructure() : BaseInfrastructure<>("test_Kronecker_creation.yml") {};
};

struct VertexState {
    double clustering_global = -1;
    double diameter = -1;
    double mean_deg = -1;
    double var = -1;
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

// -- Tests -------------------------------------------------------------------
// Undirected graphs
BOOST_FIXTURE_TEST_CASE(create_Kron_graph, Test_Graph)
{
    test_config_callable (

      [&](auto test_cfg){

        const auto g0 = GraphCreation::create_Kronecker_graph<Graph>(test_cfg, *rng, logger);
        std::size_t num_vertices = 1;
        double mean_degree = 1;
        bool is_complete = true;

        for (const auto& factor_map : test_cfg["Kronecker"]) {

            const auto factor_cfg = factor_map.second;
            const auto model = get_as<std::string>("model", factor_cfg);
            const auto N = get_as<std::size_t>("num_vertices", factor_cfg);

            num_vertices *= N;
            if (model == "complete") {
                mean_degree *= N;
            }
            else if (model == "chain") {
                double k = 2*(N-1)/(static_cast<double>(N));
                mean_degree *= (k+1);
                is_complete = false;
            }
            else {
                mean_degree *= (get_as<double>("mean_degree", factor_cfg)+1);
                is_complete = false;
            }
        }

        mean_degree -= 1;
        const auto num_edges = std::round(static_cast<double>(num_vertices)/2*mean_degree);

        BOOST_TEST(boost::num_vertices(g0) == num_vertices);
        BOOST_TEST(boost::num_edges(g0) == num_edges);

        TestUtils::assert_no_parallel_self_edges(g0);

        if (is_complete) {
            for (const auto& v : Utopia::range<Utopia::IterateOver::vertices>(g0)) {
                double c = boost::clustering_coefficient(g0, v);
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
            const auto g = AuxGraphs::create_zero_c_graph<Graph>(N, k);
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

            TestUtils::assert_no_parallel_self_edges(g);
        }
    }
}
