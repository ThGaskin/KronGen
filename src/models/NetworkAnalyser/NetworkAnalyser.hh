#ifndef UTOPIA_MODELS_NETWORKANALYSER_HH
#define UTOPIA_MODELS_NETWORKANALYSER_HH

// standard library includes
#include <random>
#include <iterator>

// third-party library includes
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/clustering_coefficient.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/range.hpp>

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/graph.hh>
#include <utopia/data_io/graph_utils.hh>

#include "graph_metrics.hh"

namespace Utopia::Models::NetworkAnalyser
{
// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
using vector = typename std::vector<double>;

/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<>;

// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename GraphType>
class NetworkAnalyser : public Model<NetworkAnalyser<GraphType>, ModelTypes>
{
  public:
    /// The type of the Model base class of this derived class
    using Base = Model<NetworkAnalyser, ModelTypes>;

    /// Data type of the group to write model data to, holding datasets
    using DataGroup = typename Base::DataGroup;

    /// Data type for a dataset
    using DataSet = typename Base::DataSet;

    // .. Graph-related types and rule types ..................................
    /// Data type for a vertex descriptor
    using VertexDesc =
        typename boost::graph_traits<GraphType>::vertex_descriptor;

    using vertices_size_type = typename boost::graph_traits<GraphType>::vertices_size_type;


    /// Data type for an edge descriptor
    using EdgeDesc = typename boost::graph_traits<GraphType>::edge_descriptor;

  private:

    // -- Members -------------------------------------------------------------
    /// A re-usable uniform real distribution to evaluate probabilities
    std::uniform_real_distribution<double> _prob_distr;

    /// The graph
    GraphType _g;

    /// The diameter of the graph
    std::size_t _diam;

    // .. Datagroups ..........................................................

    // Which features to write
    const std::pair<std::string, bool> _betweenness;
    const std::pair<std::string, bool> _closeness;
    const std::pair<std::string, bool> _clustering_coeff;
    const std::pair<std::string, bool> _core_number;
    const std::pair<std::string, bool> _degree;
    const std::pair<std::string, bool> _diameter;
    const std::pair<std::string, bool> _distance_avg;
    const std::pair<std::string, bool> _distance_harmonic;
    const std::pair<std::string, bool> _distance_max;
    const std::pair<std::string, bool> _reciprocity;

    /// Graph datagroup
    const std::shared_ptr<DataGroup> _dgrp_g;

    // .. Datasets ............................................................
    const std::shared_ptr<DataSet> _dset_betweenness;
    const std::shared_ptr<DataSet> _dset_closeness;
    const std::shared_ptr<DataSet> _dset_clustering_coeff;
    const std::shared_ptr<DataSet> _dset_core_number;
    const std::shared_ptr<DataSet> _dset_degree;
    const std::shared_ptr<DataSet> _dset_diameter;
    const std::shared_ptr<DataSet> _dset_distance_avg;
    const std::shared_ptr<DataSet> _dset_distance_harmonic;
    const std::shared_ptr<DataSet> _dset_distance_max;
    const std::shared_ptr<DataSet> _dset_reciprocity;

  public:
    // -- Model Setup ---------------------------------------------------------

    template<class ParentModel>
    NetworkAnalyser (
        const std::string& name,
        ParentModel& parent_model,
        GraphType& graph,
        const DataIO::Config& custom_cfg = {}
    )
    :
        // Initialize first via base model
        Base(name, parent_model, custom_cfg),

        // Initialize the uniform real distribution to range [0., 1.]
        _prob_distr(0., 1.),

        // Now initialize the graph
        _g {graph},

        _betweenness("betweenness", get_as<bool>("centralities",
                                                this->_cfg["graph_analysis"])),
        _closeness("closeness", get_as<bool>("centralities",
                                                this->_cfg["graph_analysis"])),
        _clustering_coeff("clustering_coeff", get_as<bool>("clustering_coeff",
                                                this->_cfg["graph_analysis"])),
        _core_number("core_number", get_as<bool>("core_number",
                                                this->_cfg["graph_analysis"])),
        _degree("degree", get_as<bool>("degree", this->_cfg["graph_analysis"])),
        _diameter("diameter", get_as<bool>("diameter", this->_cfg["graph_analysis"])),
        _distance_avg("distance_avg", get_as<bool>("distances",
                                                this->_cfg["graph_analysis"])),
        _distance_harmonic("distance_harmonic", get_as<bool>("distances",
                                                this->_cfg["graph_analysis"])),
        _distance_max("distance_max", get_as<bool>("distances",
                                                this->_cfg["graph_analysis"])),
        _reciprocity("reciprocity", get_as<bool>("reciprocity",
                                                this->_cfg["graph_analysis"])),

        _dgrp_g(create_graph_group(_g, this->_hdfgrp, "graph_data")),

        // Datasets
        _dset_betweenness(this->create_dataset(_betweenness)),
        _dset_closeness(this->create_dataset(_closeness)),
        _dset_clustering_coeff(this->create_dataset(_clustering_coeff)),
        _dset_core_number(this->create_dataset(_core_number)),
        _dset_degree(this->create_dataset(_degree)),
        _dset_diameter(this->create_dataset(_diameter)),
        _dset_distance_avg(this->create_dataset(_distance_avg)),
        _dset_distance_harmonic(this->create_dataset(_distance_harmonic)),
        _dset_distance_max(this->create_dataset(_distance_max)),
        _dset_reciprocity(this->create_dataset(_reciprocity))

    {
        if (_betweenness.second) {
            _dset_betweenness->add_attribute("is_vertex_property", true);
            _dset_betweenness->add_attribute("dim_name__1", "vertex_idx");
            _dset_betweenness->add_attribute("coords_mode__vertex_idx", "trivial");
        }
        if (_closeness.second) {
            _dset_closeness->add_attribute("is_vertex_property", true);
            _dset_closeness->add_attribute("dim_name__1", "vertex_idx");
            _dset_closeness->add_attribute("coords_mode__vertex_idx", "trivial");
        }
        if (_clustering_coeff.second) {
            _dset_clustering_coeff->add_attribute("is_vertex_property", true);
            _dset_clustering_coeff->add_attribute("dim_name__1", "vertex_idx");
            _dset_clustering_coeff->add_attribute("coords_mode__vertex_idx", "trivial");
        }
        if (_core_number.second) {
            _dset_core_number->add_attribute("is_vertex_property", true);
            _dset_core_number->add_attribute("dim_name__1", "vertex_idx");
            _dset_core_number->add_attribute("coords_mode__vertex_idx", "trivial");
        }
        if (_degree.second) {
            _dset_degree->add_attribute("dim_name__1", "vertex_idx");
            _dset_degree->add_attribute("coords_mode__vertex_idx", "trivial");
        }
        // if (_diameter.second) {
        //     _dset_diameter->add_attribute("is_vertex_property", true);
        //     _dset_diameter->add_attribute("dim_name__0", "time");
        // }
        if (_distance_avg.second) {
            _dset_distance_avg->add_attribute("dim_name__1", "vertex_idx");
            _dset_distance_avg->add_attribute("coords_mode__vertex_idx", "trivial");
            _dset_distance_harmonic->add_attribute("dim_name__1", "vertex_idx");
            _dset_distance_harmonic->add_attribute("coords_mode__vertex_idx", "trivial");
            _dset_distance_max->add_attribute("dim_name__1", "vertex_idx");
            _dset_distance_max->add_attribute("coords_mode__vertex_idx", "trivial");
        }
        if (_reciprocity.second) {
            _dset_reciprocity->add_attribute("dim_name__1", "vertex_idx");
            _dset_reciprocity->add_attribute("coords_mode__vertex_idx", "trivial");
        }

        this->_log->info("Saving the graph ...");

        save_graph(_g, _dgrp_g);

        this->_log->info("Graph saved.");
  }

  private:

    // Only initialize datasets for entries selected in the cfg
    std::shared_ptr<DataSet> create_dataset(
        const std::pair<std::string, bool>& selection)
    {
        if (get_as<bool>("enabled", this->_cfg["graph_analysis"])
            and selection.second)
        {
            if (selection.first == "diameter") {
              return this->create_dset(
                  selection.first, _dgrp_g, {}
              );
            }
            else {
                return this->create_dset(
                    selection.first, _dgrp_g, {boost::num_vertices(_g)}
                );
            }
        }
        else {
            return 0;
        }
    }


  public:

    void prolog() {

      if (get_as<bool>("enabled", this->_cfg["graph_analysis"])) {

          const auto num_vertices = boost::num_vertices(_g);

          // Ensure connectedness
          if (_distance_avg.second) {

              this->_log->info("Checking for connectedness ... ");
              size_t counter = 0;
              for (const auto v : range<IterateOver::vertices>(_g)) {
                  if (boost::degree(v, _g) == 0) {
                     auto w = random_vertex(_g, *this->_rng);
                     while (w == v) {
                       w = random_vertex(_g, *this->_rng);
                     }
                     boost::add_edge(v, w, _g);
                     ++counter;
                  }
              }
              this->_log->info("Done: added {} edges.", counter);
          }

          this->_log->info("Starting graph analysis ...");

          vector betweenness(num_vertices);
          vector closeness(num_vertices);
          size_t deg = 0;
          std::vector<VertexDesc> vertices(num_vertices);


          std::vector<std::vector<size_t>> D(num_vertices);

          if (_betweenness.second) {
              this->_log->info("Calculating centralities ...");
              std::pair<vector, vector> centralities =
                  get_centralities(_g);
              betweenness = centralities.first;
              closeness = centralities.second;
              this->_log->info("Done.");
          }

          for (const auto v : range<IterateOver::vertices>(_g)) {

              if (v%(num_vertices/100)==0){
                  this->_log->info("Iterating ... status: {}\%", 100*static_cast<double>(v)/num_vertices);
              }

              deg = boost::degree(v, _g);
              D[deg].push_back(v);

              if (_betweenness.second) {
                _g[v].state.betweenness = betweenness[v];
                _g[v].state.closeness = closeness[v];
              }

              if (_clustering_coeff.second) {
                _g[v].state.clustering_coeff = clustering_coefficient(_g, v);
              }

              if (_degree.second) {
                _g[v].state.degree = deg;
              }

              if (_distance_avg.second) {
                  get_distances(_g, v, num_vertices);
              }

          }

          if (_core_number.second) {
              this->_log->info("Computing core numbers ... ");

              compute_core_numbers(_g, D);

          }

          if (_diameter.second) {
              this->_log->info("Computing the diameter ... ");
              auto starting_point = fourSweep<vertices_size_type>(_g);
               _diam = iFUB(starting_point.first, starting_point.second, 0, _g);
          }

          this->_log->info("Graph analysis complete.");

        }
    }
    void perform_step() {}
    void monitor() {}
    void write_data() {}
    void epilog() {

      this->_log->info ("Writing data ... ");

      auto [v, v_end] = boost::vertices(_g);

      if (_betweenness.second) {
          _dset_betweenness->write(v, v_end, [this](const auto v) {
              return this->_g[v].state.betweenness;
          });
          _dset_closeness->write(v, v_end, [this](const auto v) {
              return this->_g[v].state.closeness;
          });
      }

      if (_clustering_coeff.second){
          _dset_clustering_coeff->write(v, v_end, [this](const auto v) {
              return this->_g[v].state.clustering_coeff;
          });
      }

      if (_core_number.second){
          _dset_core_number->write(v, v_end, [this](const auto v) {
              return this->_g[v].state.core_number;
          });
      }

      if (_degree.second){
          _dset_degree->write(v, v_end, [this](const auto v) {
              return this->_g[v].state.degree;
          });
      }

      if (_diameter.second){
          _dset_diameter->write(_diam);
      }

      if (_distance_avg.second){
          _dset_distance_avg->write(v, v_end, [this](const auto v) {
              return this->_g[v].state.distance_avg;
          });
          _dset_distance_harmonic->write(v, v_end, [this](const auto v) {
              return this->_g[v].state.distance_harmonic;
          });
          _dset_distance_max->write(v, v_end, [this](const auto v) {
              return this->_g[v].state.distance_max;
          });
      }

      this->_log->info ("Data written.");

    }

};

}  // namespace Utopia::Models::NetworkAnalyser

#endif  // UTOPIA_MODELS_NETWORKANALYSER_HH
