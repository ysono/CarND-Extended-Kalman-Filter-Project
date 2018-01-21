#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class KalmanFilter {
public:

  // state vector
  VectorXd x_;

  // state covariance matrix
  MatrixXd P_;

  // identity matrix with the dimension of `x_`'s
  MatrixXd I_;

  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param x Initial state
   */
  void Init(const VectorXd &x);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict(const MatrixXd &F,
               const MatrixXd &Q);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const VectorXd &z,
              const MatrixXd &H,
              const MatrixXd &R);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const VectorXd &z,
                 const VectorXd &z_pred,
                 const MatrixXd &Hj,
                 const MatrixXd &R);
};

#endif /* KALMAN_FILTER_H_ */
