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
#include "../graph_properties.hh"
#include "../graph_types.hh"
#include "../KronGen.hh"
#include "../utils.hh"
#include "test_utils.hh"

using namespace Utopia::TestTools;
using namespace Utopia::Models::KronGen;


// -- Types -------------------------------------------------------------------
// Initialise the core logger
auto logger = Utopia::init_logger("root.KronGen", spdlog::level::info);


struct Infrastructure : public BaseInfrastructure<> {
    Infrastructure() : BaseInfrastructure<>("test_Kronecker_creation.yml") {};
};

/// The test graph types
struct Test_Graph : Infrastructure {

  // undirected
  using Graph = Utopia::Models::KronGen::GraphType;

};

// -- Tests -------------------------------------------------------------------
// Undirected graphs
BOOST_FIXTURE_TEST_CASE(create_Kron_graph, Test_Graph)
{
    test_config_callable (

      [&](auto test_cfg){

        const auto g0 = GraphCreation::create_graph<Graph>(test_cfg, *rng, logger, false);
        std::size_t num_vertices = 1;
        double mean_degree = 0;
        int tensor_power = -1;
        bool is_complete = true;
        try {
            tensor_power = get_as<int>("power", test_cfg["Kronecker"]);
        }
        catch (YAML::InvalidNode&){}
        catch (Utopia::KeyError&){}

        if (tensor_power > 0) {
            const auto& factor_map = get_as<Config>("Kronecker", test_cfg);
            const auto& factor_cfg = get_as<Config>("Graph1", factor_map);
            const auto N = get_as<std::size_t>("num_vertices", factor_cfg);
            const auto model = get_as<std::string>("model", factor_cfg);
            const auto t = GraphTypes::to_graphtype(model);
            double k = 0;
            try {
                k = get_as<double>("mean_degree", factor_cfg);
            }
            catch (YAML::InvalidNode&){}
            catch (Utopia::KeyError&){}
            k = GraphProperties::mean_degree(N, k, t);

            is_complete = (model == "complete");

            num_vertices = pow(N, tensor_power);
            mean_degree = pow(k+1, tensor_power)-1;
        }
        else {
            for (const auto& factor_map : test_cfg["Kronecker"]) {

                const auto factor_cfg = factor_map.second;
                const auto model = get_as<std::string>("model", factor_cfg);
                is_complete = (is_complete and (model == "complete"));
                const auto N = get_as<std::size_t>("num_vertices", factor_cfg);
                const auto t = GraphTypes::to_graphtype(model);
                double k = 0;
                try {
                    k = get_as<double>("mean_degree", factor_cfg);
                }
                catch (YAML::InvalidNode&){}
                catch (Utopia::KeyError&){}
                k = GraphProperties::mean_degree(N, k, t);

                num_vertices = Utils::Kronecker_num_vertices(N, num_vertices);
                mean_degree = Utils::Kronecker_mean_degree(mean_degree, k);
            }
        }

        const auto num_edges = std::round(static_cast<double>(num_vertices)/2*mean_degree);

        BOOST_TEST(boost::num_vertices(g0) == num_vertices);
        BOOST_TEST(boost::num_edges(g0) == num_edges);

        TestUtils::assert_no_parallel_self_edges(g0);

        if (is_complete) {
            BOOST_TEST(num_vertices == mean_degree+1);
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
