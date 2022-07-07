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

std::vector<double> sosfilt(const std::vector<double>& sos, const std::vector<double>& x) {
    ///TODO:

}

/**
 * Низкочастотный. На входе data, частота среза и частота дискретизации
 */
std::vector<double> butter_filter(const std::vector<double>& data, const double fs, const double ffs, auto order = 4) {
    ///TODO:
    const auto wn = 2 * fs / ffs;
    std::vector<double> result;
//    sos = butter(order, wn, btype='low', output='sos')
//    filtered_values = sosfilt(sos, data)
    return result;
}

/**
 * Полосовой. На входе data, частота дискретизации, lowcut и highcut
 */
std::vector<double> butter_bandpass(const std::vector<double>& data, const double fs, const double lowcut, const double highcut, auto order = 4) {
    ///TODO:
    const auto nyq = 0.5 * fs;
    const auto low = lowcut / nyq;
    const auto high = highcut / nyq;
//    sos = butter(order, [low, high], btype='band', output='sos')
    std::vector<double> filtered;
//    filtered_values = sosfilt(sos, data)
    return filtered;
}

#endif //SAFETY_CALCULATION_DATA_PROCESSING_HPP
