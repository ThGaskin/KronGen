#ifndef UTOPIA_MODELS_KRONGEN_GRAPHTYPES
#define UTOPIA_MODELS_KRONGEN_GRAPHTYPES

namespace Utopia::Models::KronGen::GraphTypes{

// All possible graph types considered
enum GraphType {
    Chain,
    Complete,
    ErdosRenyi,
    KlemmEguiluz,
    Regular,
    ScaleFree,
    SmallWorld
};

// Convenient printing function for graph types
std::string Graph_Type[] = {
    "Chain",
    "Complete",
    "Erdos-Renyi",
    "Regular",
    "Scale-free",
    "Small-world"
};

// Convert string to GraphType
GraphType to_graphtype(std::string s) {
    if (s == "BarabasiAlbert" or s == "scale-free"){
        return ScaleFree;
    }
    else if (s == "chain") {
        return Chain;
    }
    else if (s == "complete") {
        return Complete;
    }
    else if (s == "ErdosRenyi") {
        return ErdosRenyi;
    }
    else if (s == "KlemmEguiluz") {
        return KlemmEguiluz;
    }
    else if (s == "regular" or s == "Regular") {
        return Regular;
    }
    else if (s == "WattsStrogatz"){
        return SmallWorld;
    }
    else {
        return ErdosRenyi;
    }
}

} // namespace KronGen::GraphTypes

#endif // UTOPIA_MODELS_KRONGEN_GRAPHTYPES
