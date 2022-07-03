#ifndef SAFETY_CALCULATION_GPS_HPP
#define SAFETY_CALCULATION_GPS_HPP

#include <cmath>
#include <vector>
#include <deque>
#include "data.hpp"


/**
 * Формула гаверсинуса
 */
double get_distance_between_points(double lat1, double lon1, double lat2, double lon2) {
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
 * Медиана массива
 */
double median(auto data) {
    std::sort(data.begin(), data.end());
    if (data.size() % 2 != 0) {
        return data[data.size() / 2];
    }
    return (data[(data.size() - 1) / 2] + data[data.size() / 2]) / 2;
}

/**
 * Медианный фильтр по одной оси
 */
std::vector<double> median_filter(const std::vector<double>& data, const std::size_t window) {
    if (data.size() < window) {
        return {};
    }

    std::vector<double> result;
    result.reserve(data.size() -  (window - 1) + 2);

    std::deque<double> frame;
    frame.emplace_back(data[0]);
    for (auto iter = data.begin(); iter < data.begin() + window - 1; ++iter) {
        frame.emplace_back(*iter);
    }
    for (std::size_t i = window - 1; i < data.size(); ++i) {
        result.emplace_back(median(frame));
        frame.pop_front();
        frame.emplace_back(data[i]);
    }
    result.emplace_back(median(frame));
    frame.pop_front();
    frame.emplace_back(data[data.size() - 1]);
    result.emplace_back(median(frame));

    return result;
}

void calculate_distance(GPS& gps) {
    const std::size_t length = gps.lat.size();
    std::vector<double> distance(length);
    double prev_lat = 0;
    double prev_lon = 0;
    for (std::size_t i = 0; i < length; ++i) {
        if (prev_lon == 0) {
            distance.emplace_back(0);
        } else {
            distance.emplace_back(get_distance_between_points(gps.lat[i],gps.lon[i],prev_lat, prev_lon));
        }
        prev_lon = gps.lon[i];
        prev_lat = gps.lat[i];
    }
    gps.filtered = median_filter(distance, 3);
    gps.route.reserve(gps.filtered.size());
    double route = 0;
    for (std::size_t i = 0; i < gps.filtered.size(); ++i) {
        route += gps.filtered[i];
        gps.route.emplace_back(route);
    }
}

std::vector<std::array<double, 4>> get_equations(std::vector<double>& time,  std::vector<double>& lat, const double acc_time) {
    std::vector<std::array<double,4>> result;
    if (acc_time < time[0]) {
        time.emplace(time.begin(), acc_time);
        lat.emplace(time.begin(), 0);
    }
    for (std::size_t i = 1; i < time.size(); ++i) {
        if (time[i - 1] != time[i]) {
            const double a = (lat[i - 1] - lat[i]) / (time[i - 1] - time[i]);
            const double b = lat[i - 1] - a * time[i - 1];
            result.emplace_back(std::array{time[i - 1], time[i], a, b});
        }
    }
    return result;
}

std::vector<double> get_x_dist(std::vector<double>& time, std::vector<double>& lat, const std::vector<double>& acc_time) {
    auto equations = get_equations(time, lat, acc_time[0]);
    std::vector<double> y;
    std::size_t index = 0;
    for (const auto& time : acc_time) {
        if (time > equations[index][1]) {
            index += 1;
        }
        if (index < equations.size()) {
            y.emplace_back(equations[index][2] * time + equations[index][3]);
        }
        if (index == equations.size()) {
            break;
        }
    }

    return y;
}

#endif //SAFETY_CALCULATION_GPS_HPP
