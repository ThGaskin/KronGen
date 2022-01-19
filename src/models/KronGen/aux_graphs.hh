#ifndef UTOPIA_MODELS_KRONGEN_AUX_GRAPHS
#define UTOPIA_MODELS_KRONGEN_AUX_GRAPHS

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/random.hpp>
#include <boost/property_map/dynamic_property_map.hpp>

#include "utopia/core/graph/iterator.hh"
#include "utopia/core/types.hh"
#include "utopia/core/graph/creation.hh"

#include "graph_types.hh"

namespace Utopia::Models::KronGen::AuxGraphs {

/// Auxiliary graphs used in the KronGen model

/// Create a k-regular graph with zero clustering.
/**
  * \param N    Number of vertices. Must be at least 2*k, and even.
  * \param k    Mean degree of graph. Must be at least 2.
  * \return g   The graph
*/
template<typename Graph>
Graph create_zero_c_graph(const std::size_t N, const std::size_t k)
{

    if (N < 2*k) {
        throw std::invalid_argument("N must be greater or equal to 2k!");
    }
    if (N%2) {
        throw std::invalid_argument("N must be even!");
    }
    if (k < 2) {
        throw std::invalid_argument("k must be greater or equal to 2!");
    }

    Graph g{N};

    for (std::size_t v = 0; v < N/2; ++v) {
        const auto v_2 = v + N/2 - k;
        for (std::size_t i = 0; i < k; ++i) {
            auto nb = (v_2 + i)%(N/2) + N/2;
            add_edge(v, nb, g);
        }
    }

    return g;
}

/// Create a 1-d lattice chain of length N. The diameter of the chain is N-1
/**
  * \param N   Length of chain
*/
template<typename Graph>
Graph create_chain_graph(const std::size_t N)
{
    if (N == 2) {
        return Utopia::Graph::create_regular_graph<Graph>(N, 1, false);
    }
    Graph g = Utopia::Graph::create_regular_graph<Graph>(N, 2, false);
    remove_edge(0, 1, g);

    return g;
}
/// Create a star-graph with a given number of vertices, a given mean_degree,
/// and a given diameter
/**
  * \param N          Number of vertices
  * \param k          mean degree
  * \param diameter   Diameter
  *
  * \returns g        Graph
*/
// To do: Randomise this!
template<typename Graph, typename RNGType>
Graph create_star_graph(const std::size_t N,
                        const double k,
                        const std::size_t diameter,
                        RNGType& rng)
{
    if ((N*k-2*diameter+2) > (N-diameter+1)*(N-diameter)) {
        throw std::invalid_argument("Star graph cannot be created!");
    }

    std::uniform_real_distribution<double> distr(0, 1);
    Graph g{N};
    const std::size_t center = 0;
    for (std::size_t i = 1; i < N; ++i) {
        add_edge(center, i, g);
    }

    // At this point, the graph has diameter 2 and mean degree 2*(N-2)/N
    auto target = k * N;
    while (2*num_edges(g) < target) {
        auto w = random_vertex(g, rng);
        auto v = random_vertex(g, rng);
        while ((w == v) or (edge(v, w, g).second)){
            w = random_vertex(g, rng);
            v = random_vertex(g, rng);
        }
        add_edge(w, v, g);
    }

    // The graph now has diameter 2 and mean degree k
    if (diameter > 2) {
        // Create the anchor path of length diam
        std::vector<std::size_t> anchor_path;
        auto start = random_vertex(g, rng);
        while (start == center) {
            start = random_vertex(g, rng);
        }
        anchor_path.push_back(start);
        std::size_t i = 0;
        std::size_t edges_to_add = 0;

        while (i < (diameter - 2)) {

            if (i == 0) {
                edges_to_add += (degree(start, g)-1);
            }
            else if (degree(start, g) > 2) {
                edges_to_add += (degree(start, g)-2);
            }
            auto v = random_vertex(g, rng);
            while ((v == center)
                or (v == start)
                or (std::find(anchor_path.begin(), anchor_path.end(), v) != anchor_path.end()))
            {
                v = random_vertex(g, rng);
            }
            clear_vertex(start, g);
            add_edge(start, v, g);
            if (anchor_path.back() != start) {
              add_edge(start, anchor_path.back(), g);
            }

            if (i != 0) {
                anchor_path.push_back(start);
            }

            start = v;

            ++i;
        }

        for (std::size_t j = 0; j < edges_to_add; ++j) {
            auto x = random_vertex(g, rng);
            auto y = random_vertex(g, rng);
            while ((x == y)
                or (std::find(anchor_path.begin(), anchor_path.end(), x) != anchor_path.end())
                or (std::find(anchor_path.begin(), anchor_path.end(), y) != anchor_path.end())
                or (edge(x, y, g).second)){
                    x = random_vertex(g, rng);
                    y = random_vertex(g, rng);
            }
            add_edge(x, y, g);
        }
    }

    return g;

}

template<typename Graph, typename RNGType>
Graph create_graph(const size_t N,
                   const size_t k,
                   const GraphTypes::GraphType t,
                   RNGType& rng,
                   double c_t) {

    if ((N == k+1) or (t == GraphTypes::Complete)) {
        return Utopia::Graph::create_complete_graph<Graph>(N);
    }
    else if (t == GraphTypes::Chain) {
        return create_chain_graph<Graph>(N);
    }
    else if (t == GraphTypes::ErdosRenyi) {
        return Utopia::Graph::create_ErdosRenyi_graph<Graph>(N, k, false, false, rng);
    }
    else if (t == GraphTypes::KlemmEguiluz) {
        return Utopia::Graph::create_KlemmEguiluz_graph<Graph>(N, k, 1-c_t, rng);
    }
    else if (t == GraphTypes::Regular) {
        return Utopia::Graph::create_regular_graph<Graph>(N, k, false);
    }
    else {
        return Graph{};
    }
}
/// Extended create_graph function: calls the Utopia::Graph function of the
/// same name
template<typename Graph, typename RNGType>
Graph create_graph(const Config& graph_cfg, RNGType& rng)
{
    const std::string model = get_as<std::string>("model", graph_cfg);

    if (model == "zero_c") {
        return create_zero_c_graph<Graph>(
            get_as<std::size_t>("num_vertices", graph_cfg),
            get_as<std::size_t>("mean_degree", graph_cfg)
        );
    }

    else if (model == "chain") {
        return create_chain_graph<Graph>(
            get_as<std::size_t>("num_vertices", graph_cfg)
        );
    }

    else if (model == "star") {
        return create_star_graph<Graph>(
            get_as<std::size_t>("num_vertices", graph_cfg),
            get_as<double>("mean_degree", graph_cfg),
            get_as<std::size_t>("diameter", graph_cfg["star"]),
            rng
        );
    }

    else {
        return Utopia::Graph::create_graph<Graph>(graph_cfg, rng);
    }
}

}
#endif
