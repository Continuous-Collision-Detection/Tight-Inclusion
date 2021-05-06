#pragma once

#include <vector>

#include <Eigen/Core>

namespace inclusion_ccd {

Eigen::MatrixXd
read_rational_csv(const std::string& inputFileName, std::vector<bool>& results);

} // namespace ccd
