#ifndef UTOPIA_MODELS_KRONGEN_TESTUTILS
#define UTOPIA_MODELS_KRONGEN_TESTUTILS

#include <boost/mpl/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

namespace Utopia::Models::KronGen::TestUtils {

// -- Helper function ----------------------------------------------------------
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

}
#endif
