#ifndef SAFETY_CALCULATION_DATA_HPP
#define SAFETY_CALCULATION_DATA_HPP

#include <vector>

struct Accel {
    std::vector<double> cart_t;
    std::vector<double> cart_v;
    std::vector<double> cab_y_t;
    std::vector<double> cab_y_v;
    std::vector<double> cab_z_t;
    std::vector<double> cab_z_v;
};

struct GPS {
    std::vector<double> time;
    std::vector<double> lat;
    std::vector<double> lon;
    std::vector<double> speed;
};

#endif //SAFETY_CALCULATION_DATA_HPP
