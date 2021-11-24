#ifndef UTOPIA_MODELS_NETWORKANALYSER__GRAPH_METRICS_HH
#define UTOPIA_MODELS_NETWORKANALYSER__GRAPH_METRICS_HH

// standard library includes
#include <random>
#include <iterator>

// third-party library includes
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/clustering_coefficient.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/properties.hpp>
#include <boost/range.hpp>

#include "diameter.hh"

namespace Utopia::Models::NetworkAnalyser
{

using vector = typename std::vector<double>;

/// Calculate the betweenness centrality of each vertex.
/// Calculate the closeness/relative betweenness centrality for each vertex
/// (normalized with the highest possible value which would be reached
///  if a node is crossed by every single shortest path).
template<typename GraphType>
const std::pair<vector, vector> get_centralities(const GraphType& g)
{
    vector centrality(num_vertices(g));

    boost::brandes_betweenness_centrality(g,
          boost::make_iterator_property_map(centrality.begin(),
                                       get(boost::vertex_index, g),
                                       double())
    );

    vector closeness = centrality;

    boost::relative_betweenness_centrality(g,
            boost::make_iterator_property_map(closeness.begin(),
                                              get(boost::vertex_index, g),
                                              double())
    );

    return std::make_pair(centrality, closeness);
}

/// Calculate the global clustering coefficient
template <typename GraphType>
double global_clustering_coeff(const GraphType& g)
{
    if (num_vertices(g) == 2 and num_edges(g) == 1) {
      return 0;
    }
    std::size_t num_triangles = 0;
    std::size_t degree_sum = 0;
    for (const auto v : range<IterateOver::vertices>(g)) {
        // Undirected: factor of 2
        num_triangles += 2*boost::num_triangles_on_vertex(g, v);
        degree_sum += boost::degree(v, g)*(degree(v, g)-1);
    }

    return static_cast<double>(num_triangles)/degree_sum;

}

// Identify groups of agents that are connected via out-edges.
// Note that completely isolated vertices are also identified
// as closed community.
template<typename GraphType>
const std::vector<std::vector<size_t>> closed_communities(const GraphType& g) {

    std::vector<std::vector<size_t>> cc;
    std::vector<size_t> temp_c;
    bool next;

    // Find all communities through a loop over all vertices as the
    // source of the community.
    for (const auto v : range<IterateOver::vertices>(g)) {

        next = false;

        // If vertex is part of an already discovered community
        // its community has to be the same (if out-degree > 0).
        for (auto& c: cc) {
            if (std::find(c.begin(), c.end(), v) != c.end()) {
                next = true;
            }
        }

        if (next) {
            continue;
        }

        else {
            if (boost::in_degree(v, g) < 2) {
                // This is the case of a 'loner'.
                for (auto& c: cc) {
                    for (auto& w: c) {
                        if (boost::edge(v, w, g).second) {
                            c.push_back(v);
                            next = true;
                            break;
                        }
                    }
                }
            }
            else {
                // Else get the community originating from the vertex.
                temp_c.clear();
                fill_community(v, temp_c, g);
                cc.push_back(temp_c);
            }
        }
    }

    return cc;
}


// Computes the k-core number of each vertex. The core number of a vertex is the
// lowest number number n such that when vertices with degree <= n are removed
// iteratively (that is, including those that have degree <= n once other
// vertices have already been removed), that vertex remains in the graph.
// This runs in linear time.
template<typename GraphType>
void compute_core_numbers(GraphType& g,
                          std::vector<std::vector<size_t>> D)
{

    GraphType h;
    boost::copy_graph(g, h);

    size_t N = num_vertices(h);
    size_t k = 0;
    for (size_t n = 0; n < N; ++n) {
      size_t i = 0;
      while (D[i].empty()) { ++i; }
      const auto v = D[i].back();
      k = std::max(k, i);
      g[v].state.core_number = k;

      for (const auto e : range<IterateOver::out_edges>(v, h)) {
        auto w = target(e, h);
        size_t curr_deg = boost::in_degree(w, h);
          D[curr_deg].erase(std::remove(D[curr_deg].begin(), D[curr_deg].end(), w), D[curr_deg].end());
          D[curr_deg-1].push_back(w);
      }
      boost::clear_vertex(v, h);
      D[i].erase(std::find(D[i].begin(), D[i].end(), v));
    }
}


/// Get various indicators of shortest path lengths: average,
/// harmonic average, max
template<typename GraphType, typename VertexDesc>
void get_distances(
    GraphType& g,
    VertexDesc v,
    const std::size_t num_vertices)
{
    // BFS of all vertex distances
    auto d = get_distances(v, g);

    const double max = *std::max_element(d.begin(), d.end());
    const double sum = std::accumulate(d.begin(), d.end(), 0.0);
    const double avg = sum / (num_vertices - 1);

    double harmonic = 0.;
    for (auto && elem: d) {
        elem != 0 ? harmonic += 1. / elem : harmonic += 0;
    }

    g[v].state.distance_avg = avg;
    harmonic != 0
        ? g[v].state.distance_harmonic = 1 / (harmonic / (num_vertices -1))
        : g[v].state.distance_harmonic = 0;
    g[v].state.distance_max = max;

}

// Calculate the reciprocity for a single node (= fraction of outgoing
// links for which the mutual link exists as well).
// For directed graphs only.
template<typename Graph, typename VertexDescType>
double reciprocity(const VertexDescType v, const Graph& g)
{
    double r = 0.;
    for (const auto w : range<IterateOver::neighbors>(v, g)) {
        if (edge(w, v, g).second) {
            ++r;
        }
    }

    return r / out_degree(v, g);
}

// Calculate the reciprocity of the whole graph (= fraction of mutual links).
// For directed graphs only.
template<typename Graph>
double reciprocity(const Graph& g)
{
    double r = 0;
    for (const auto e : range<IterateOver::edges>(g)) {
        if (edge(target(e, g), source(e, g), g).second)
        {
            ++r;
        }
    }

    return r / num_edges(g);
}

/// Calculate the mean degree and the degree variance of a graph
// To do: test this
template <typename GraphType>
std::pair<double, double> degree_statistics (const GraphType& g) {
    double e = 0;
    double ee = 0;
    for (const auto v : range<IterateOver::vertices>(g)) {
        const auto d = degree(v, g);
        e += d;
        ee += pow(d, 2);
    }

    const std::size_t N = num_vertices(g);
    const double k = e/N;
    const double var = (ee/N - pow(e/N, 2));

    return {k, var};
}

} // namespace Utopia::Models::NetworkAnalyser


#endif // UTOPIA_MODELS_NETWORKANALYSER__GRAPH_METRICS_HH
