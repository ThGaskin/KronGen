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
#include "../grid_search.hh"
#include "../type_definitions.hh"
#include "../utils.hh"

using namespace Utopia;
using namespace Utopia::Models::KronGen::GridSearch;
using namespace Utopia::Models::KronGen::Utils;
using namespace Utopia::Models::KronGen::TypeDefinitions;
using namespace Utopia::TestTools;

struct Infrastructure : public BaseInfrastructure<> {
    Infrastructure() : BaseInfrastructure<>("test_grid_search.yml") {};
};

/// The test graph types
struct Test_Graph : Infrastructure {

  // undirected
  using Graph = Utopia::Models::KronGen::NWType;

};

// Test getting simple factors
BOOST_AUTO_TEST_CASE (test_get_N_factors)
{
    const std::vector<size_t> targets = {1, 2, 17, 52, 169, 500, 5000};
    for (const auto& N : targets) {
        auto factors = get_factors(N, false);
        for (const auto& f : factors) {
            BOOST_TEST(f.size() == 2);
            BOOST_TEST(f[0] * f[1] == N);
            BOOST_TEST(f[0] != 1);
            BOOST_TEST(f[1] != 1);
        }
        factors = get_factors(N, true);
        BOOST_TEST(factors.size() != 0);
        BOOST_TEST((factors[0] == std::vector<size_t>{1, N}));
        for (const auto& f : factors) {
            BOOST_TEST(f.size() == 2);
            BOOST_TEST(f[0] * f[1] == N);
        }
    }

    BOOST_CHECK_THROW (get_factors(0, true), std::invalid_argument);
    BOOST_CHECK_THROW (get_factors(0, false), std::invalid_argument);

}

BOOST_AUTO_TEST_CASE (test_get_k_factors)
{
    const std::vector<size_t> targets = {1, 2, 6, 16, 50, 51, 169, 500, 5000};
    for (const auto& k : targets) {
        auto factors = get_k_factors(k, false);
        for (const auto& f : factors) {
            BOOST_TEST(f.size() == 2);
            BOOST_TEST((f[0]+1) * (f[1]+1) == (k+1));
            BOOST_TEST(f[0] > 0);
            BOOST_TEST(f[1] > 0);
        }
        factors = get_k_factors(k, true);
        BOOST_TEST(factors.size() != 0);
        BOOST_TEST((factors[0] == std::vector<size_t>{0, k}));
        for (const auto& f : factors) {
            BOOST_TEST(f.size() == 2);
            BOOST_TEST((f[0]+1) * (f[1]+1) == (k+1));
        }
    }
}

BOOST_AUTO_TEST_CASE (test_get_d_factors)
{
      const std::vector<size_t> N = {2, 17, 52, 169, 500, 5000};
      for (const auto& n : N){
            for (size_t d = 1; d < 10; ++d){
                auto fac = get_d_factors(n, get_factors, d);
                for (const auto& f : fac) {
                    BOOST_TEST(f.size() == d);
                    size_t res = 1;
                    for (const auto& ff: f) {
                        res *= ff;
                        if (f.size() > 1){
                            BOOST_TEST(ff > 1);
                        }
                    }
                    BOOST_TEST(res == n);
                }
            }
      }

      for (const auto& k : N){
            for (size_t d = 0; d < 10; ++d){
                auto fac = get_d_factors(k, get_k_factors, 10);
                for (const auto& f : fac) {
                    BOOST_TEST(f.size() == d);
                    size_t res = 1;
                    for (const auto& ff: f) {
                        res *= (ff+1);
                        BOOST_TEST(ff >= 0);
                    }
                    BOOST_TEST(res == k+1);
                }
            }
      }

      // Test an explicit case
      auto test = get_d_factors (52, get_factors, 1);
      factors res = {{52}};
      BOOST_TEST(test == res);
      test = get_d_factors (52, get_factors, 2);
      res = {{2, 26}, {4, 13}};
      BOOST_TEST(test == res);
      test = get_d_factors (52, get_factors, 3);
      res = {{2, 2, 13}};
      BOOST_TEST(test == res);
}


// Test getting the grid in N
BOOST_AUTO_TEST_CASE (test_N_grid)
{
    const std::vector<size_t> targets = {2, 17, 52, 169, 500, 5000};
    const std::vector<double> errors = {0., 0.1, 0.2};
    const std::vector<size_t> min_dim = {1, 1, 1, 1, 2, 3, 6};
    const std::vector<double> max_dim = {1, 2, 3, 4, 2, 5, 6};
    for (const auto& N_target : targets){
        for (const auto& err : errors) {
            for (size_t d = 0; d < min_dim.size(); ++d){
                const auto N_grid = get_N_grid (N_target, err, min_dim[d], max_dim[d]);
                for (const auto& n : N_grid) {
                    size_t N_res = 1;
                    for (const auto& fac : n) {
                        BOOST_TEST(fac != 1);
                        N_res *= fac;
                    }
                    BOOST_TEST(1.0*N_res >= (N_target*(1.-err)));
                    BOOST_TEST(1.0*N_res <= (N_target*(1.+err)));
                    BOOST_TEST(n.size() >= min_dim[d]);
                    BOOST_TEST(n.size() <= max_dim[d]);
                }
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

// Test the validity check
BOOST_AUTO_TEST_CASE (test_validity_check)
{
    map_type test_map{{"degree_sequence",
      entry_type{{"calculate", {true, vector_pt{{0, 1}}}},
                 {"target", {false, ""}}}
    }};
    std::vector<size_t> N = {10, 12, 17, 19, 31};
    std::vector<size_t> k = {2, 4, 8, 5, 2};
    BOOST_TEST(check_validity(N, k, test_map));

    test_map["degree_sequence"]["target"].first = true;
    BOOST_TEST(!check_validity(N, k, test_map));

    k[0] = 3;
    k[4] = 3;
    BOOST_TEST(check_validity(N, k, test_map));

    k.emplace_back(4);
    BOOST_TEST(!check_validity(N, k, test_map));

    N.emplace_back(3);
    BOOST_TEST(!check_validity(N, k, test_map));

    N[5] = 5;
    BOOST_TEST(check_validity(N, k, test_map));

    k[5] = 1;
    BOOST_TEST(!check_validity(N, k, test_map));

    N[5] = 2;
    BOOST_TEST(check_validity(N, k, test_map));

    k[5] = 0;
    N[5] = 1;
    BOOST_TEST(check_validity(N, k, test_map));

    N = {};
    k = {};
    BOOST_TEST(!check_validity(N, k, test_map));

}

// Test finding possible graph types
BOOST_AUTO_TEST_CASE (test_type_search)
{
    map_type test_map = map_type{
      {"degree_sequence", entry_type{{"calculate", {true, vector_pt{{0, 1}}}},
        {"target", {false, ""}}}}
    };
    test_map.insert({"diameter", entry_type{{"target", {true, 30.0}}}});

    std::vector<size_t> N = {10, 12, 17, 19, 31};
    std::vector<size_t> k = {2, 4, 8, 5, 2};

    auto types = find_possible_types(N, k, test_map, 0.1);

    BOOST_TEST((types[0] == std::vector<GraphType>{N.size(), GraphType::ErdosRenyi}));
    bool found_chain = false;
    for (const auto& type : types) {
        if (type.back() == GraphType::Chain) {
            found_chain = true;
        }
    }
    BOOST_TEST(found_chain);
}

// // // Test the grid search
// BOOST_FIXTURE_TEST_CASE (test_grid_search, Test_Graph)
// {
//
// }
