#define BOOST_TEST_MODULE config extraction test

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <spdlog/spdlog.h>

#include <utopia/core/testtools.hh>
#include <utopia/core/types.hh>
#include <utopia/core/graph.hh>

#include "../KronGen.hh"
#include "../graph_creation.hh"
#include "../utils.hh"

using namespace Utopia;
using namespace Utopia::Models::KronGen::GraphCreation;
using namespace Utopia::Models::KronGen::GraphProperties;
using namespace Utopia::Models::KronGen::Utils;
using namespace Utopia::TestTools;

struct Infrastructure : public BaseInfrastructure<> {
    Infrastructure() : BaseInfrastructure<>("test_config_extraction.yml") {};
};

/// The test graph types
struct Test_Graph : Infrastructure {

  // undirected
  using Graph = Utopia::Models::KronGen::NWType;

};

// Test extraction of analysis targets from config
BOOST_FIXTURE_TEST_CASE(test_Kronecker_extraction_1, Test_Graph)
{
    const auto config = get_as<Config>("params", cfg["Kronecker1"]);
    const auto analysis_cfg = get_as<Config>("graph_analysis", config["NetworkAnalyser"]);
    const auto analysis_targets = get_analysis_targets(analysis_cfg);

    BOOST_TEST(analysis_targets.size() == 6);
    for (auto& p : analysis_targets) {
        BOOST_TEST(
            (p.first == "clustering"
             or p.first == "degree_sequence"
             or p.first == "degree_variance"
             or p.first == "diameter"
             or p.first == "mean_degree"
             or p.first == "num_vertices")
        );
    }
}

BOOST_FIXTURE_TEST_CASE(test_Kronecker_extraction_2, Test_Graph)
{
    const auto config = get_as<Config>("params", cfg["Kronecker2"]);
    const auto analysis_cfg = get_as<Config>("graph_analysis", config["NetworkAnalyser"]);
    const auto analysis_targets = get_analysis_targets(analysis_cfg);

    BOOST_TEST(analysis_targets.size() == 3);
    for (auto& p : analysis_targets) {
        BOOST_TEST(
            (p.first == "degree_sequence"
             or p.first == "diameter"
             or p.first == "num_vertices")
        );
    }
}

BOOST_FIXTURE_TEST_CASE(test_KronGen_extraction, Test_Graph)
{
    const auto config = get_as<Config>("params", cfg["KronGen"]);
    const auto analysis_cfg = get_as<Config>("graph_analysis", config["NetworkAnalyser"]);
    auto targets_cfg = get_as<Config>("KronGen", config["create_graph"]);
    targets_cfg = get_as<Config>("targets", targets_cfg);
    const auto analysis_targets = get_analysis_targets(analysis_cfg, targets_cfg);

    BOOST_TEST(analysis_targets.size() == 12);
}

// Test calculation of properties
BOOST_FIXTURE_TEST_CASE(test_Kronecker_calculation, Test_Graph) {
    const auto config = get_as<Config>("params", cfg["Kronecker1"]);
    const auto analysis_cfg = get_as<Config>("graph_analysis", config["NetworkAnalyser"]);
    auto analysis_targets = get_analysis_targets(analysis_cfg);

    const auto G = Utopia::Graph::create_complete_graph<Graph>(20);

    calculate_properties(G, analysis_targets, log);

    BOOST_TEST(std::any_cast<double>(
        analysis_targets["num_vertices"]["calculate"].second) == 20
    );
    BOOST_TEST(std::any_cast<double>(
        analysis_targets["mean_degree"]["calculate"].second) == 19
    );
    BOOST_TEST(std::any_cast<double>(
        analysis_targets["diameter"]["calculate"].second) == 1
    );
    BOOST_TEST(std::any_cast<double>(
        analysis_targets["clustering"]["calculate"].second) == 1
    );
    BOOST_TEST(std::any_cast<double>(
        analysis_targets["degree_variance"]["calculate"].second) == 0
    );
    BOOST_TEST((std::any_cast<vector_pt>(
        analysis_targets["degree_sequence"]["calculate"].second) == vector_pt{{19, 20}})
    );

}

// Test properties are written
BOOST_FIXTURE_TEST_CASE(test_KronGen_write, Test_Graph) {
    const auto config = get_as<Config>("params", cfg["Write"]);
    const auto analysis_cfg = get_as<Config>("graph_analysis", config["NetworkAnalyser"]);
    auto analysis_targets = get_analysis_targets(analysis_cfg);

    auto G = Utopia::Graph::create_complete_graph<Graph>(100);

    calculate_properties(G, analysis_targets, log);

    write_properties(G, analysis_targets);

    BOOST_TEST(G[0].state.num_vertices == 100);
    BOOST_TEST(G[0].state.mean_degree == 99);
    BOOST_TEST(G[0].state.clustering_global == 1);
    BOOST_TEST(G[0].state.diameter == 1);
    BOOST_TEST((G[0].state.degree_sequence == vector_pt{{99, 100}}));
    BOOST_TEST(G[0].state.degree_variance == 0);

}
