#ifndef SAFETY_CALCULATION_GPS_HPP
#define SAFETY_CALCULATION_GPS_HPP

#include <cmath>

/**
 * Формула гаверсинуса
 */
double get_distance_between_points(double lat1, double lon1, double lat2, double lon2) {
    ///TODO: ensure parameters are in radians
    const double d_lon = lon2 - lon1;
    const double d_lat = lat2 - lat1;
    const double p
            = std::pow(std::sin(d_lat / 2), 2.0)
            + std::cos(lat1) * std::cos(lat2) * std::pow(std::sin(d_lon / 2), 2.0);
    const double q = 2 * std::asin(std::sqrt(p));
    const double EARTH_RADIUS = 6371.0;
    return q * EARTH_RADIUS;
}

/**
 * Медианный фильтр по одной оси
 */
double median_filter(/*data, n*/) {
    ///TODO:
}

double calculate_distance(/*data*/) {
    ///TODO:
}

///*?*/ get_equations(/*data, time*/) {
//    ///TODO:
//}
//
///*?*/ get_x_dist(/*gps_data, acc_time*/) {
//    ///TODO:
//}

#endif //SAFETY_CALCULATION_GPS_HPP
