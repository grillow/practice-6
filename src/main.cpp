#include <iostream>
#include "DataReader.hpp"
#include <memory>

void print_acc_info(const Accel& acc) {
    std::cout << "Данные от акселерометров:\n";
    std::cout << "Время КП:\t" << acc.cart_t.size() << std::endl;
    std::cout << "Данные КП:\t" << acc.cart_v.size() << std::endl;
    std::cout << "Время Y:\t" << acc.cab_y_t.size() << std::endl;
    std::cout << "Данные Y:\t" << acc.cab_y_v.size() << std::endl;
    std::cout << "Время Z:\t" << acc.cab_z_t.size() << std::endl;
    std::cout << "Данные Z:\t" << acc.cab_z_v.size() << std::endl;
}

auto main(int argc, char **argv) -> int {
    auto reader = std::make_unique<CSVDataReader>("data/accel.csv", "data/gps.csv");

    auto acc = reader->readAccel();
    print_acc_info(acc);

    auto gps = reader->readGPS();

    return 0;
}
