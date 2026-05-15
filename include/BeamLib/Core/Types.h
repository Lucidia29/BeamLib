#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>

namespace beamlib {

using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Mat2 = Eigen::Matrix2d;
using Mat3 = Eigen::Matrix3d;
using VecX = Eigen::VectorXd;
using MatX = Eigen::MatrixXd;
using SpMat = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

template <int N>
using VecN = Eigen::Matrix<double, N, 1>;

template <int M, int N>
using MatMN = Eigen::Matrix<double, M, N>;

using Tensor4_3x3 = std::array<std::array<Mat3, 3>, 3>;

} // namespace beamlib
