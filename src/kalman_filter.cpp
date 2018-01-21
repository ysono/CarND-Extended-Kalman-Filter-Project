#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(const VectorXd &x) {
  x_ = x;
  long x_size = x_.size();
  I_ = MatrixXd::Identity(x_size, x_size);
}

void KalmanFilter::Predict(const MatrixXd &F,
                           const MatrixXd &Q) {
  x_ = F * x_;
  P_ = F * P_ * F.transpose() + Q;
}

void KalmanFilter::Update(const VectorXd &z,
                          const MatrixXd &H,
                          const MatrixXd &R) {
  VectorXd h = H * x_;
  UpdateEKF(z, h, H, R);
}

void KalmanFilter::UpdateEKF(const VectorXd &z,
                             const VectorXd &h,
                             const MatrixXd &Hj,
                             const MatrixXd &R) {
  VectorXd y = z - h;
  MatrixXd Hjt = Hj.transpose();
  MatrixXd S = Hj * P_ * Hjt + R;
  MatrixXd K = P_ * Hjt * S.inverse();

  x_ = x_ + K * y;
  P_ = (I_ - K * Hj) * P_;
}
