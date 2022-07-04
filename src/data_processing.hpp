#ifndef SAFETY_CALCULATION_DATA_PROCESSING_HPP
#define SAFETY_CALCULATION_DATA_PROCESSING_HPP

#include "data.hpp"
#include <vector>
#include <cmath>

std::vector<double> get_zones(const std::vector<double>& route, const double length) {
    double lim = length;
    double zone = 0;
    std::vector<double> zones;

    for (auto& r : route) {
        if (r > lim) {
            lim += length;
            ++zone;
        }
        zones.emplace_back(zone);
    }

    return zones;
}



#endif //SAFETY_CALCULATION_DATA_PROCESSING_HPP
