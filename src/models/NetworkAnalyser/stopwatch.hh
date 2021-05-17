#ifndef UTOPIA_MODELS_NETWORKANALYSER_HH_STOPWATCH_HH
#define UTOPIA_MODELS_NETWORKANALYSER_HH_STOPWATCH_HH



// standard library includes
#include <chrono>
#include <utility> // tuple
#include <string>
#include <map>
#include <vector>


// third-party library includes



// Utopia-related includes



// ContDiseaseNet includes
//#include "helpers.hh"



namespace Utopia::Models::NetworkAnalyser {


using name = std::string;
using time_point
        = std::chrono::time_point<std::chrono::high_resolution_clock>;
using duration = std::chrono::duration<double>;
using values = std::tuple<double, time_point, time_point>;
using Timer = std::map<name, values>;
using Name = std::map<name, name>;


class Stopwatch {
public:

private:
    Timer timers;
    Name names;

public:
    void start(const name& key, const name& full = "");
    void stop(const name& key);
    void start_next(const name& key, const name& full = "");
    std::vector<std::string> get_keys() const;
    std::vector<name> get_names() const;
    std::vector<double> get_durations() const;
private:

};


void Stopwatch::start(const name& key, const name& full) {
    auto full_name = (full == "") ? key : full;
    std::cout<<"Running "<<full<<"..."<<std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    auto val = std::make_tuple(0., start, start);
    if ( timers.find(key) != timers.end() ) {
        std::get<1>(timers[key]) = start;
    } else {
        timers.insert(std::make_pair(key, val));
        names.insert(std::make_pair(key, full_name));
    }
}

void Stopwatch::stop(const name& key) {
    auto stop = std::chrono::high_resolution_clock::now();
    duration dur = stop - std::get<1>(timers[key]);
    std::get<0>(timers[key]) += dur.count();
    std::get<2>(timers[key]) = stop;
}

void Stopwatch::start_next(const name& key, const name& full) {
    auto now = std::chrono::high_resolution_clock::now();

    if (std::get<0>(timers.rbegin()->second) == 0.) {
        std::cout << "passed this\n\n";
        duration dur = now - std::get<1>(timers.rbegin()->second);
        std::get<0>(timers.rbegin()->second) += dur.count();
        std::get<2>(timers.rbegin()->second) = now;
    }
    auto full_name = (full == "") ? key : full;
    auto val = std::make_tuple(0., now, now);
    timers.insert(std::make_pair(key, val));
    names.insert(std::make_pair(key, full_name));
}

std::vector<std::string> Stopwatch::get_keys() const {
    std::vector<std::string> keys_vector;
    for (Name::const_iterator it = names.begin();
                                it != names.end(); ++it) {
        keys_vector.push_back(it->first);
    }
    return keys_vector;
}

std::vector<name> Stopwatch::get_names() const {
    std::vector<name> names_vector;
    for (Name::const_iterator it = names.begin();
                                it != names.end(); ++it) {
        names_vector.push_back(it->second);
    }
    return names_vector;
}

std::vector<double> Stopwatch::get_durations() const {
    std::vector<double> durations_vector;
    for (Timer::const_iterator it = timers.begin();
                                it != timers.end(); ++it) {
        durations_vector.push_back(std::get<0>(it->second));
    }
    return durations_vector;
}




Stopwatch stopw;





} // namespace Utopia::Models::ContDiseaseNet


#endif // UTOPIA_MODELS_CONTDISEASENET_STOPWATCH_HH
