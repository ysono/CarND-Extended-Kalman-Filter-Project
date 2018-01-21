#ifndef TOOLS_H_
#define TOOLS_H_
#include <tuple>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  VectorXd CalcRadarX(const VectorXd &z);
  VectorXd CalcLaserX(const VectorXd &z);

  std::tuple<VectorXd, MatrixXd> CalcRadarH(const VectorXd &x_state);
  void CorrectRadarHPhase(const VectorXd &z, VectorXd &h);

  VectorXd CalculateRMSE(VectorXd &previous_sse,
                         int &previous_sse_count,
                         const VectorXd &estimation,
                         const VectorXd &ground_truth);
};

#endif /* TOOLS_H_ */
