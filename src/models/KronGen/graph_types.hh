#ifndef UTOPIA_MODELS_KRONGEN_GRAPHTYPES
#define UTOPIA_MODELS_KRONGEN_GRAPHTYPES

namespace Utopia::Models::KronGen::GraphTypes{

enum GraphType {
    Chain,
    Complete,
    ErdosRenyi,
    Regular
};

std::string Graph_Type[] = {
    "Chain",
    "Complete",
    "ErdosRenyi",
    "Regular"
};

}
#endif
