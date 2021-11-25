#ifndef UTOPIA_MODELS_KRONGEN_GRAPHTYPES
#define UTOPIA_MODELS_KRONGEN_GRAPHTYPES

namespace Utopia::Models::KronGen::GraphTypes{

// All possible graph types considered
enum GraphType {
    Chain,
    Complete,
    ErdosRenyi,
    Regular
};

// Convenient printing function for graph types
std::string Graph_Type[] = {
    "Chain",
    "Complete",
    "Erdos-Renyi",
    "Regular"
};

} // namespace KronGen::GraphTypes

#endif // UTOPIA_MODELS_KRONGEN_GRAPHTYPES
