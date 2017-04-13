#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include <stdexcept>
#include "Eigen/Dense"

class Tools {
private:
  Tools() {};

public:
  /**
  * A helper method to calculate RMSE.
  */
  static Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  static Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

  static double normalize_angle(double angle);
};

#define ZERO (1e-6)
#define PI 3.1415926535897932384626433832795

#endif /* TOOLS_H_ */
