#ifndef SAFETY_CALCULATION_DATA_PROCESSING_HPP
#define SAFETY_CALCULATION_DATA_PROCESSING_HPP

#include "data.hpp"
#include <vector>
#include <cmath>
#include <complex>
#include <numbers>

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

// https://github.com/scipy/scipy/blob/4cf21e753cf937d1c6c2d2a0e372fbc1dbbeea81/scipy/signal/_filter_design.py#L4256
struct buttap_result {
    std::vector<double> z;
    std::vector<std::complex<double>> p;
    double k;
};
buttap_result buttap(const std::size_t N) {
    using namespace std::complex_literals;

    std::vector<std::complex<double>> p;
    p.reserve(N);
    for (int i = 0, value = static_cast<int>(-N) + 1; i < N; ++i) {
        const auto power = 1i * std::numbers::pi * static_cast<double>(value)
                / (2.0 * static_cast<double>(N));
        const auto exp = -std::exp(power);
        p.emplace_back(exp);
        value += 2;
    }

    return {
        .z = {},
        .p = std::move(p),
        .k = 1
    };
}

void butter_low_sos(const std::size_t order, const double wn) {
    ///TODO:
    // buttap

}

/**
 * Низкочастотный. На входе data, частота среза и частота дискретизации
 */
std::vector<double> butter_filter(const std::vector<double>& data, const double fs, const double ffs, const std::size_t order = 4) {
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
