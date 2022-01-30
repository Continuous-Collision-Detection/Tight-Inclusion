#pragma once

#include <vector>

#include <Eigen/Core>

namespace ticcd::rational {

    Eigen::MatrixXd read_rational_csv(
        const std::string &inputFileName, std::vector<bool> &results);
    Eigen::MatrixXd
    read_csv(const std::string &inputFileName, std::vector<double> &results);

} // namespace ticcd::rational
