#ifndef SAFETY_CALCULATION_DATA_PROCESSING_HPP
#define SAFETY_CALCULATION_DATA_PROCESSING_HPP

#include "data.hpp"
#include <vector>
#include <cmath>

double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}

void gps_deg2rad(GPS& gps) {
    for (auto &e : gps.lat) e = deg2rad(e);
    for (auto &e : gps.lon) e = deg2rad(e);
}

#endif //SAFETY_CALCULATION_DATA_PROCESSING_HPP
