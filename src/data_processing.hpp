#ifndef SAFETY_CALCULATION_DATA_PROCESSING_HPP
#define SAFETY_CALCULATION_DATA_PROCESSING_HPP

#include "data.hpp"
#include <vector>
#include <cmath>
#include <complex>
#include <numbers>
#include <numeric>

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

// https://github.com/scipy/scipy/blob/4cf21e753cf937d1c6c2d2a0e372fbc1dbbeea81/scipy/signal/_filter_design.py
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

ZPK lp2lp_zpk(ZPK zpk, const double wo = 1.0) {
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

ZPK lp2bp_zpk(ZPK zpk, const double wo = 1.0, const double bw = 1.0) {
    for (auto& z : zpk.z) {
        z *= bw / 2;
    }
    for (auto& p : zpk.p) {
        p *= bw / 2;
    }

    zpk.z.reserve(zpk.z.size() * 2);
    const auto z_initial_size = zpk.z.size();
    for (std::size_t i = 0; i < z_initial_size; ++i) {
        const auto z = zpk.z[i];
        const auto first = z + std::sqrt(z * z - wo * wo);
        const auto second = z - std::sqrt(z * z - wo * wo);
        zpk.z[i] = first;
        zpk.z.emplace_back(second);
    }
    zpk.p.reserve(zpk.p.size() * 2);
    const auto p_initial_size = zpk.p.size();
    for (std::size_t i = 0; i < p_initial_size; ++i) {
        const auto p = zpk.p[i];
        const auto first = p + std::sqrt(p * p - wo * wo);
        const auto second = p - std::sqrt(p * p - wo * wo);
        zpk.p[i] = first;
        zpk.p.emplace_back(second);
    }

    const auto degree = zpk.p.size() - zpk.z.size();
    zpk.z.reserve(zpk.z.size() + degree);
    for (std::size_t i = 0; i < degree; ++i) {
        zpk.z.emplace_back(0.0);
    }

    zpk.k *= std::pow(bw, degree);

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

// matrix: N rows, 6 cols
using SOS = std::vector<std::array<double, 6>>;

SOS zpk2sos(ZPK zpk) {
    ///TODO:
    // pairing = 'nearest'
    if (zpk.z.size() == zpk.p.size() == 0) {
        ///TODO:
//        return np.array([[k, 0., 0., 1., 0., 0.]])
    }

    const auto zp_max_size = std::max(zpk.z.size(), zpk.p.size());
    zpk.z.reserve(zp_max_size);
    zpk.p.reserve(zp_max_size);
    for (std::size_t i = zpk.z.size(); i < zp_max_size; ++i) {
        zpk.z.emplace_back(0);
    }
    for (std::size_t i = zpk.p.size(); i < zp_max_size; ++i) {
        zpk.p.emplace_back(0);
    }

    const std::size_t n_sections = (zp_max_size + 1) / 2;

    // possible optimization: reserve more size in the previous step
    if (zpk.p.size() % 2 == 1) {
        zpk.p.emplace_back(0);
        zpk.z.emplace_back(0);
    }

    std::sort(zpk.z.begin(), zpk.z.end());

    // might be incorrect:
    // 1. remove numbers with negative imaginary part
    // 2. sort first by real part, and then by magnitude of imaginary part
    zpk.p.erase(std::remove_if(zpk.p.begin(), zpk.p.end(),
                               [](const std::complex<double>& e){ return e.imag() < 0; }
                               ), zpk.p.end());
    std::sort(zpk.p.begin(), zpk.p.end(), [](const auto& left, const auto& right){
        if (left.real() != right.real()) {
            return left.real() < right.real();
        }
        return left.imag() < right.imag();
    });

    const auto idx_worst = [](const std::vector<std::complex<double>>& p) {
        return std::distance(p.cbegin(),
         std::min_element(p.cbegin(), p.cend(), [](const auto& left, const auto& right) {
                return
                    std::abs(1 - std::abs(left))
                    <
                    std::abs(1 - std::abs(right));
            })
         );
    };

    SOS sos(n_sections); // zero initialized

    std::vector<std::complex<double>> z;
    z.reserve(zpk.z.size());
    for (const auto& e : zpk.z) {
        z.emplace_back(e);
    }

    // [n_sections - 1 .. 0]
    for (auto si = n_sections - 1; si >= 0; --si) {
        const auto single_zpk2sos = [](const auto& z, const auto& p, const auto& k) {
            const auto poly = [](const std::vector<std::complex<double>>& vec) {
                ///TODO:
            };

            std::vector<double> sos(6);

            ///TODO:
            // b, a = zpk2tf(z, p, k)
            // sos = [b, a]

            // b = k * poly(z)
            // a = poly(p)

            return sos;
        };

        const auto nearest_real_idx
                = [](const std::vector<std::complex<double>>& fro, const std::complex<double>& to) -> std::size_t {
            std::vector<std::size_t> order;
            order.reserve(fro.size());
            for (std::size_t i = 0; i < fro.size(); ++i) {
                order.emplace_back(i);
            }
            std::sort(order.begin(), order.end(), [&fro, &to](std::size_t left, std::size_t right) {
                return std::abs(fro[left] - to)
                <
                std::abs(fro[right] - to);
            });

            std::vector<bool> mask;
            mask.reserve(fro.size());
            for (const auto& e : fro) {
                mask.emplace_back(e.imag() == 0);
            }

            std::vector<std::size_t> nonzero;
            for (const auto& e : mask) {
                nonzero.emplace_back(e != 0);
            }

            return order[nonzero[0]];
        };
        const auto nearest_complex_idx
                = [](const std::vector<std::complex<double>>& fro, const std::complex<double>& to) -> std::size_t {
                    std::vector<std::size_t> order;
                    order.reserve(fro.size());
                    for (std::size_t i = 0; i < fro.size(); ++i) {
                        order.emplace_back(i);
                    }
                    std::sort(order.begin(), order.end(), [&fro, &to](std::size_t left, std::size_t right) {
                        return std::abs(fro[left] - to)
                               <
                               std::abs(fro[right] - to);
                    });

                    std::vector<bool> mask;
                    mask.reserve(fro.size());
                    for (const auto& e : fro) {
                        mask.emplace_back(e.imag() != 0);
                    }

                    std::vector<std::size_t> nonzero;
                    for (const auto& e : mask) {
                        nonzero.emplace_back(e != 0);
                    }

                    return order[nonzero[0]];
                };
        const auto nearest_any_idx
                = [](const std::vector<std::complex<double>>& fro, const std::complex<double>& to) -> std::size_t {
                    std::vector<std::size_t> order;
                    order.reserve(fro.size());
                    for (std::size_t i = 0; i < fro.size(); ++i) {
                        order.emplace_back(i);
                    }
                    std::sort(order.begin(), order.end(), [&fro, &to](std::size_t left, std::size_t right) {
                        return std::abs(fro[left] - to)
                               <
                               std::abs(fro[right] - to);
                    });

                    return order[0];
                };

        const auto p1_idx = idx_worst(zpk.p);
        const auto p1 = zpk.p[p1_idx];
        zpk.p.erase(zpk.p.begin() + p1_idx);

        const auto p_isreal_sum = std::accumulate(zpk.p.cbegin(), zpk.p.cend(), 0, [](std::size_t sum, std::complex<double> a){
            return sum + (a.imag() == 0);
        });
        const auto z_isreal_sum = std::accumulate(z.cbegin(), z.cend(), 0, [](std::size_t sum, std::complex<double> a){
            return sum + (a.imag() == 0);
        });
        if (p1.imag() == 0 && p_isreal_sum == 0) {
            const auto z1_idx = nearest_real_idx(z, p1);
            const auto z1 = z[z1_idx];
            z.erase(z.begin() + z1_idx);
            ///TODO:
//            sos[si] = _single_zpksos([z1, 0], [p1, 0], 1)
        } else if (zpk.p.size() + 1 == z.size() && p1.imag() != 0 && p_isreal_sum == 1 && z_isreal_sum == 1) {
            const auto z1_idx = nearest_complex_idx(z, p1);
            const auto z1 = z[z1_idx];
            z.erase(z.begin() + z1_idx);
            ///TODO:
//            sos[si] = _single_zpksos([z1, z1.conj()], [p1, p1.conj()], 1)
        } else {
            std::size_t p2_idx;
            std::complex<double> p2;
            if (p1.imag() == 0) {
                std::vector<std::size_t> prealidx;
                std::vector<std::complex<double>> preal;
                for (std::size_t i = 0; i < zpk.p.size(); ++i){
                    if (zpk.p[i].imag() == 0) {
                        prealidx.emplace_back(i);
                        preal.emplace_back(zpk.p[i]);
                    }
                }

                p2_idx = prealidx[idx_worst(preal)];
                p2 = zpk.p[p2_idx];
                zpk.p.erase(zpk.p.begin() + p2_idx);
            } else {
                p2 = conj(p1);
            }

            if (!z.empty()) {
                const auto z1_idx = nearest_any_idx(z, p1);
                const auto z1 = z[z1_idx];
                z.erase(z.begin() + z1_idx);
                if (z1.imag() != 0) {
                    ///TODO:
//                    sos[si] = _single_zpksos([z1, z1.conj()], [p1, p2], 1)
                } else {
                    if (!z.empty()) {
                        const auto z2_idx = nearest_real_idx(z, p1);
                        const auto z2 = z[z2_idx];
                        z.erase(z.begin() + z2_idx);
                        ///TODO:
//                        sos[si] = _single_zpksos([z1, z2], [p1, p2], 1)
                    } else {
                        ///TODO:
//                        sos[si] = _single_zpksos([z1], [p1, p2], 1)
                    }
                }

            } else {
                ///TODO:
//                sos[si] = _single_zpksos([], [p1, p2], 1)
            }
        }

    }
    ///TODO:
//# put gain in first sos
//    sos[0][:3] *= k
    return sos;
}

// wn=wn, btype='low', output='sos'
SOS butter_low(const std::size_t order, const double wn) {
    auto zpk = buttap(order);
    const auto fs = 2.0;
    const auto warped = 2 * fs * std::tan(std::numbers::pi * wn / fs);
    zpk = lp2lp_zpk(zpk, warped);
    zpk = bilinear_zpk(zpk, fs);

    return zpk2sos(zpk);
}

// wn=[low, high], btype='band', output='sos'
SOS butter_band(const std::size_t order, const std::array<double, 2> wn) {
    auto zpk = buttap(order);
    const auto fs = 2.0;
    const std::array<double, 2> warped = {
        2 * fs * std::tan(std::numbers::pi * wn[0] / fs),
        2 * fs * std::tan(std::numbers::pi * wn[1] / fs)
    };
    const auto bw = warped[1] - warped[0];
    const auto wo = std::sqrt(warped[0] * warped[1]);
    zpk = lp2bp_zpk(zpk, wo, bw);
    zpk = bilinear_zpk(zpk, fs);

    return zpk2sos(zpk);
}

/*output*/void sosfilter(/*args*/) {
    ///TODO:

}

/**
 * Низкочастотный. На входе data, частота среза и частота дискретизации
 */
std::vector<double> butter_filter(const std::vector<double>& data, const double fs, const double ffs, const std::size_t order = 4) {
    ///TODO:
    const auto wn = 2 * fs / ffs;
    const auto sos = butter_low(order, wn);
//    filtered = sosfilter(sos, data)
//    return filtered;
}

/**
 * Полосовой. На входе data, частота дискретизации, lowcut и highcut
 */
std::vector<double> butter_bandpass(const std::vector<double>& data, const double fs, const double lowcut, const double highcut, const std::size_t order = 4) {
    ///TODO:
    const auto nyq = 0.5 * fs;
    const auto low = lowcut / nyq;
    const auto high = highcut / nyq;
    const auto sos = butter_band(order, {low, high});
//    filtered = sosfilter(sos, data)
//    return filtered;
}

#endif //SAFETY_CALCULATION_DATA_PROCESSING_HPP
