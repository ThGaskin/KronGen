#ifndef UTOPIA_MODELS_NETWORKANALYSER__DIAMETER_HH
#define UTOPIA_MODELS_NETWORKANALYSER__DIAMETER_HH

// standard library includes
#include <random>
#include <iterator>
#include <queue>

// third-party library includes
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/range.hpp>

// Utopia graph includes
#include <utopia/core/model.hh>
#include <utopia/core/graph.hh>
#include <utopia/data_io/graph_utils.hh>

namespace Utopia::Models::NetworkAnalyser
{

// Get the distances from a vertex to all other vertices in the graph.
template<typename Graph, typename VertexSizeType>
std::vector<std::size_t> get_distances(const VertexSizeType v, const Graph& g)
{
  const std::size_t num_vertices = boost::num_vertices(g);
  std::vector<std::size_t> distances(num_vertices);
  std::vector<bool> visited(num_vertices, false);
  std::queue<std::size_t> BFSQueue;

  // Distance from a vertex to itself is 0
  visited[v] = true;
  distances[v] = 0;
  BFSQueue.push(v);

  // Iterate over neighbors, add vertex to queue if not yet visited
  while (not BFSQueue.empty()) {
    const VertexSizeType current = BFSQueue.front();
    BFSQueue.pop();

    for (const auto nb : range<IterateOver::neighbors>(current, g)) {
          if (not visited[nb]) {
              visited[nb] = true;
              distances[nb] = distances[current] + 1;
              BFSQueue.push(nb);
          }
      }
  }

  return distances;

}

// Compute the eccentricity of a vertex, which is simply its maximum distance
// to another vertex
template<typename VertexSizeType, typename Graph>
std::size_t eccentricity(const VertexSizeType v, Graph& g)
{
  const auto distances = get_distances(v, g);
  return *std::max_element(distances.begin(), distances.end());

}

// Return a vector with all vertices sorted by their distance from a given vertex
template<typename VertexSizeType, typename Graph>
std::vector<std::vector<VertexSizeType>> get_F(const VertexSizeType v,
                                               const Graph& g)
{
  const auto distances = get_distances(v, g);
  const auto ecc = *std::max_element(distances.begin(), distances.end());

  std::vector<std::vector<VertexSizeType>> F(ecc+1);
  for (const auto w : range<IterateOver::vertices>(g)) {
      F[distances[w]].emplace_back(w);
  }

  return F;

}

// Compute the maximum eccentricity of a list of vertices
template<typename VertexSizeType, typename Graph>
std::size_t max_eccentricity(const std::vector<VertexSizeType> vertex_list,
                             const Graph& g)
{
    std::size_t max_ecc = 0;
    for (const auto& v : vertex_list) {
        const auto ecc = eccentricity(v, g);
        if (ecc > max_ecc) {
            max_ecc = ecc;
        }
    }

    return max_ecc;
}

// Get the vertex with the highest degree in the graph
template<typename VertexSizeType, typename Graph>
VertexSizeType get_max_degree_vertex(Graph& g) {

  std::size_t max_deg = 0;
  VertexSizeType max_deg_vertex = 0;

  for (const auto v : range<IterateOver::vertices>(g)) {
      if (boost::degree(v, g) > max_deg) {
          max_deg = boost::degree(v, g);
          max_deg_vertex = v;
      }
  }

  return max_deg_vertex;

}


// Return the vertex furthest away from a source vertex
template<typename VertexSizeType, typename Graph>
VertexSizeType furthest_vertex_from_source(const VertexSizeType v, Graph& g) {

  const auto distances = get_distances(v, g);
  const auto max = std::max_element(distances.begin(), distances.end());

  return static_cast<VertexSizeType>(std::distance(distances.begin(), max));

}

// Get the vertex midway between two vertices v and w
template<typename VertexSizeType, typename Graph>
VertexSizeType mid_vertex(VertexSizeType v, VertexSizeType w, Graph& g) {

  const std::size_t num_vertices = boost::num_vertices(g);
  std::vector<VertexSizeType> distances(num_vertices);
  std::vector<VertexSizeType> parent_vertices(num_vertices);
  std::vector<bool> visited(num_vertices, false);
  std::queue<VertexSizeType> BFSQueue;

  visited[v] = true;
  distances[v] = 0;
  BFSQueue.push(v);

  while (not BFSQueue.empty()) {
      const VertexSizeType current = BFSQueue.front();
      BFSQueue.pop();

      for (const auto nb : range<IterateOver::neighbors>(current, g)) {
          if (not visited[nb]) {
              visited[nb] = true;
              distances[nb] = distances[current] + 1;
              parent_vertices[nb] = current;
              BFSQueue.push(nb);
          }
      }
      // Reached w: can stop
      if (visited[w]) {
          break;
      }
  }

  std::size_t mid_distance = std::round(static_cast<double>(distances[w]) / 2.);
  VertexSizeType mid_vertex = w;
  while (mid_distance--) {
      mid_vertex = parent_vertices[mid_vertex];
  }

  return mid_vertex;

}

// fourSweep algorithm: return best estimate for lower bound and starting vector
// for iFub
template<typename VertexSizeType, typename Graph>
std::pair<std::size_t, VertexSizeType> fourSweep(Graph& g) {

    const auto v_1 = get_max_degree_vertex<VertexSizeType>(g);
    const auto a_1 = furthest_vertex_from_source(v_1, g);
    const auto b_1 = furthest_vertex_from_source(a_1, g);

    const auto v_2 = mid_vertex(a_1, b_1, g);
    const auto a_2 = furthest_vertex_from_source(v_2, g);
    const auto b_2 = furthest_vertex_from_source(a_2, g);

    const auto res = mid_vertex(a_2, b_2, g);

    const auto lower_bound = std::max(eccentricity(a_1, g), eccentricity(a_2, g));

    return {res, lower_bound};
}

// The iFUB algorithm: returns the diameter of the graph
/* \param Vertex v: starting node
   \param size_t lower_bound: lower bound for the diameter
   \param size_t tolerance: accuracy for diameter estimate; for tol = 0, we
   / obtain the actual diameter
*/
template<typename VertexSizeType, typename Graph>
std::size_t iFUB(const VertexSizeType v,
                 std::size_t l,
                 const std::size_t tolerance,
                 Graph& g)
{
  const auto ecc = eccentricity(v, g);
  std::size_t lower_bound = std::max(ecc, l);
  std::size_t upper_bound = 2 * ecc;
  const auto F = get_F(v, g);

  std::size_t i = ecc;
  while (upper_bound-lower_bound > tolerance) {
      const std::size_t new_lower_bound =
          std::max(lower_bound, max_eccentricity(F[i], g));

      if (new_lower_bound > 2 * (i-1)) {
          return new_lower_bound;
      }
      else {
          lower_bound = new_lower_bound;
          upper_bound = 2 * (i-1);
      }
      i--;
  }

  return lower_bound;

}

template<typename Graph>
double diameter (Graph& G) {

    using vertices_size_type = typename boost::graph_traits<Graph>::vertices_size_type;

    const auto s = fourSweep<vertices_size_type>(G);
    const auto diam = iFUB(s.first, s.second, 0, G);

    return diam;
}

}
#endif
