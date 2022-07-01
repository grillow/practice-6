#ifndef SAFETY_CALCULATION_DATAREADER_HPP
#define SAFETY_CALCULATION_DATAREADER_HPP

#include <vector>
#include <string>
#include "data.hpp"
#include <H5Cpp.h>
#include <fstream>
#include <sstream>

class IDataReader {
public:
    virtual ~IDataReader() = default;

    virtual Accel readAccel() = 0;
    virtual GPS readGPS() = 0;
};

class HDF5DataReader : public IDataReader {
public:
    explicit HDF5DataReader(const std::string& path) : file(path.c_str(), H5F_ACC_RDONLY) {}

    Accel readAccel() override {
        ///TODO:
        // dont forget to sort and modify values like in read_acc_hdf
        throw;
    }
    GPS readGPS() override {
        ///TODO:
        // dont forget to remove zeroes like in read_gps_hdf
        throw;
    }

private:
    H5::H5File file;
};

class CSVDataReader : public IDataReader {
public:
    explicit CSVDataReader(const std::string& accel, const std::string& gps) : f_accel(accel), f_gps(gps) {}

    Accel readAccel() override {
        Accel accel;

        std::string line;
        std::getline(f_accel, line);    // ignore the very first line
        while (std::getline(f_accel, line)) {
            std::stringstream ss(line);
            std::string value;
            std::getline(ss, value, ',');
            accel.cart_t.emplace_back(std::stod(value));
            std::getline(ss, value, ',');
            accel.cart_v.emplace_back(std::stod(value));
            std::getline(ss, value, ',');
            accel.cab_y_t.emplace_back(std::stod(value));
            std::getline(ss, value, ',');
            accel.cab_y_v.emplace_back(std::stod(value));
            std::getline(ss, value, ',');
            accel.cab_z_t.emplace_back(std::stod(value));
            std::getline(ss, value, ',');
            accel.cab_z_v.emplace_back(std::stod(value));
        }

        return accel;
    }

    GPS readGPS() override {
        GPS gps;

        std::string line;
        std::getline(f_accel, line);    // ignore the very first line
        while (std::getline(f_accel, line)) {
            std::stringstream ss(line);
            std::string value;
            std::getline(ss, value, ',');
            gps.time.emplace_back(std::stod(value));
            std::getline(ss, value, ',');
            gps.lat.emplace_back(std::stod(value));
            std::getline(ss, value, ',');
            gps.lon.emplace_back(std::stod(value));
            std::getline(ss, value, ',');
            gps.speed.emplace_back(std::stod(value));
        }

        return gps;
    }

private:
    std::ifstream f_accel;
    std::ifstream f_gps;
};

#endif //SAFETY_CALCULATION_DATAREADER_HPP
