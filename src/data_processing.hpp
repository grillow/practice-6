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

/*output*/void sosfilt(/*args*/) {
    ///TODO:

}

// https://github.com/scipy/scipy/blob/4cf21e753cf937d1c6c2d2a0e372fbc1dbbeea81/scipy/signal/_filter_design.py#L4256
struct ZPK {
    std::vector<double> z;
    std::vector<std::complex<double>> p;
    double k;
};
ZPK buttap(const std::size_t N) {
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

ZPK lp2lp_zpk(ZPK zpk, const double wo=1.0) {
    for (auto& z : zpk.z) {
        z *= wo;
    }
    for (auto& p : zpk.p) {
        p *= wo;
    }
    const auto degree = zpk.p.size() - zpk.z.size();
    zpk.k *= std::pow(wo, static_cast<double>(degree));

    return zpk;
}

ZPK bilinear_zpk(ZPK zpk, const double fs) {
    const auto degree = zpk.p.size() - zpk.z.size();
    const auto fs2 = 2 * fs;
    for (auto& e : zpk.z) {
        e = (fs2 + e) / (fs2 - e);
    }
    for (auto& e : zpk.p) {
        e = (fs2 + e) / (fs2 - e);
    }
    zpk.z.reserve(zpk.z.size() + degree);
    for (std::size_t i = 0; i < degree; ++i) {
        zpk.z.emplace_back(-1.0);
    }

    double zprod = 1;
    double pprod = 1;
    for (const auto& e : zpk.z) {
        zprod *= fs2 - e;
    }
    for (const auto& e : zpk.p) {
        pprod *= fs2 - e.real();
    }
    zpk.k *= zprod / pprod;

    return zpk;
}

struct SOS {
    ///TODO:
};

SOS zpk2sos(ZPK zpk) {
    ///TODO:
}

// btype='low', output='sos'
SOS butter(const std::size_t order, const double wn) {
    auto zpk = buttap(order);
    const auto fs = 2.0;
    const auto warped = 2 * fs * std::tan(std::numbers::pi * wn / fs);
    zpk = lp2lp_zpk(zpk, warped);
    zpk = bilinear_zpk(zpk, fs);

    return zpk2sos(zpk);
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
//    std::vector<double> filtered;
//    filtered_values = sosfilt(sos, data)
//    return filtered;
}

#endif //SAFETY_CALCULATION_DATA_PROCESSING_HPP
