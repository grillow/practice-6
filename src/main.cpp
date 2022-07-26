#include <iostream>
#include <memory>
#include "DataReader.hpp"
#include "data_processing.hpp"
#include "gps.hpp"


void print_acc_info(const Accel& acc) {
    std::cout << "Данные от акселерометров:\n";
    std::cout << "Время КП:\t"  << acc.cart_t.size()    << '\n';
    std::cout << "Данные КП:\t" << acc.cart_v.size()    << '\n';
    std::cout << "Время Y:\t"   << acc.cab_y_t.size()   << '\n';
    std::cout << "Данные Y:\t"  << acc.cab_y_v.size()   << '\n';
    std::cout << "Время Z:\t"   << acc.cab_z_t.size()   << '\n';
    std::cout << "Данные Z:\t"  << acc.cab_z_v.size()   << '\n';
    std::cout << std::flush;
}

auto main(int argc, char **argv) -> int {
//    auto reader = std::make_unique<CSVDataReader>("data/csv/small/accel.csv", "data/csv/small/gps.csv");
//
//    auto acc = reader->readAccel();
//    print_acc_info(acc);
//
//    auto gps = reader->readGPS();
//    gps_deg2rad(gps);

    const std::size_t order = 4;
    const auto wn = 2 * 10 / 610.3515625;

    auto zpk = buttap(order);
    const auto fs = 2.0;
    const auto warped = 2 * fs * std::tan(std::numbers::pi * wn / fs);
    zpk = lp2lp_zpk(zpk, warped);
    zpk = bilinear_zpk(zpk, fs);

    for (const auto& e : zpk.z) {
        std::cout << e << " ";
    }
    std::cout << std::endl;
    for (const auto& e : zpk.p) {
        std::cout << e << " ";
    }
    std::cout << std::endl;
    std::cout << zpk.k << std::endl;



    return 0;
}
