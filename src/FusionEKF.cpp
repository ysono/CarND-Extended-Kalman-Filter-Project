#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;
  previous_dt_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 500, 0,
             0, 0, 0, 500;

  F_ = MatrixXd::Identity(4, 4); // valid when dt = 0

  Q_ = MatrixXd(4, 4);

  noise_ax = 9;
  noise_ay = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  VectorXd z = measurement_pack.raw_measurements_;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    VectorXd x;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      x = tools.CalcRadarX(z);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x = tools.CalcLaserX(z);
    }
    ekf_.Init(x);

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Assume single-threaded
  if (previous_dt_ != dt) {
    previous_dt_ = dt;

    float dt2 = dt * dt;
    float dt3 = dt2 * dt;
    float dt4 = dt3 * dt;

    F_(0, 2) = dt;
    F_(1, 3) = dt;

    Q_ <<  dt4 / 4 * noise_ax, 0, dt3 / 2 * noise_ax, 0,
           0, dt4 / 4 * noise_ay, 0, dt3 / 2 * noise_ay,
           dt3 / 2 * noise_ax, 0, dt2 * noise_ax, 0,
           0, dt3 / 2 * noise_ay, 0, dt2 * noise_ay;
  }

  ekf_.Predict(F_, Q_);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    VectorXd h;
    MatrixXd Hj;
    std::tie(h, Hj) = tools.CalcRadarH(ekf_.x_);
    tools.CorrectRadarHPhase(z, h);
    ekf_.UpdateEKF(z, h, Hj, R_radar_);
  } else {
    ekf_.Update(z, H_laser_, R_laser_);
  }

  // print the output
  // cout << "x_ = " << ekf_.x_ << endl;
  // cout << "P_ = " << ekf_.P_ << endl;
}
