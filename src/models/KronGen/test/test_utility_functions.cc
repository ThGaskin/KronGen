#define BOOST_TEST_MODULE Kronecker graph utils test

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/clustering_coefficient.hpp>

#include <utopia/core/types.hh>
#include <utopia/core/graph.hh>

#include "../graph_creation.hh"
#include "../../NetworkAnalyser/graph_metrics.hh"
#include "../utils.hh"

using namespace Utopia::Models::KronGen::GraphCreation;
using namespace Utopia::Models::KronGen::Utils;
using namespace Utopia::Models::NetworkAnalyser;
using namespace Utopia;

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
                std::cout<<p<<": "<<f.first<<", "<<f.second<<std::endl;
                BOOST_TEST((!((p-1) % f.first) or !((p+1) % f.first)));
                BOOST_TEST((!((p-1) % f.second) or !((p+1) % f.second)));
            }
        }

    }
}

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
    }
}

// To do: test diameter estimation, extreme cases
