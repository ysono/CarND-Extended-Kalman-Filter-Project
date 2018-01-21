#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalcRadarX(const VectorXd &z) {
  float ro = z(0);
  float phi = z(1);
  float ro_dot = z(2);

  float cos_phi = cos(phi);
  float sin_phi = sin(phi);

  VectorXd x(4);
  x << ro * cos_phi,
       ro * sin_phi,
       ro_dot * cos_phi,
       ro_dot * sin_phi;
  return x;
}

VectorXd Tools::CalcLaserX(const VectorXd &z) {
  VectorXd x(4);
  x << z(0),
       z(1),
       0,
       0;
  return x;
}

std::tuple<VectorXd, MatrixXd> Tools::CalcRadarH(const VectorXd& x_state) {
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float c1 = pow(px, 2) + pow(py, 2);
  float c2 = sqrt(c1);

  float h_rho_dot;

  MatrixXd Hj(3,4);

  if (fabs(c1) < 0.0001) {
    // If if the object is measured to be extremely close to self, then
    // - it's reasonable to apprximate rho as zero, hence set the first row of Hj to zero.
    // - it's impossible to estimate phi, so choose zero, i.e. the direction of self's travel.
    // - it's impossible to estimate rho_dot, so choose zero, i.e. no acceleration.

    h_rho_dot = 0;

    Hj << 0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0;
  } else {
    h_rho_dot = (px * vx + py * vy) / c2;

    float c3 = c1 * c2;
    float c4 = vx*py - vy*px;

    Hj << (px/c2), (py/c2), 0, 0,
          -(py/c1), (px/c1), 0, 0,
          py*c4/c3, -px*c4/c3, px/c2, py/c2;
  }

  VectorXd h(3);
  h << c2,
       atan2(py, px),
       h_rho_dot;

  return std::make_tuple(h, Hj);
}

void Tools::CorrectRadarHPhase(const VectorXd &z, VectorXd &h) {
  float z_phi = z(1);
  float h_phi = h(1);

  float diff = h_phi - z_phi;
  while (diff > M_PI) {
    h_phi -= 2 * M_PI;
    diff = h_phi - z_phi;
  }
  while (diff < -M_PI) {
    h_phi += 2 * M_PI;
    diff = h_phi - z_phi;
  }

  h(1) = h_phi;
}

VectorXd Tools::CalculateRMSE(VectorXd &previous_sse,
                              int &previous_sse_count,
                              const VectorXd &estimation,
                              const VectorXd &ground_truth) {
  auto err = (estimation - ground_truth).array();
  VectorXd sq = err * err;
  previous_sse += sq;
  previous_sse_count ++;

  VectorXd mse = previous_sse / previous_sse_count;
  VectorXd rmse = mse.array().sqrt();
  return rmse;
}

