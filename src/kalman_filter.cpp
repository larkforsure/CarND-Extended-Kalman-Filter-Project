#include <iostream>
#include "tools.h"
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    /**
     * predict the state
     */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

// for Ladar
void KalmanFilter::Update(const VectorXd &z) {
    /**
     * update the state by using Kalman Filter equations
     */
    // compare predicates with measuments
    VectorXd y = z - (H_ * x_);  // error
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    // filter out the new state values and uncertainity
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

// for Radar
void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
     * update the state by using Extended Kalman Filter equations
     */
    // compare predicates with measuments
    // have initialize H_ with Hj jacobian before calling
    VectorXd y = z - h(x_);  // error
    y(1) = Tools::normalize_angle(y(1)); // need this when use atan2
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    // filter out the new state values and uncertainity
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

// h here, Hj is caculated by Tools::CalculateJacobian 
VectorXd KalmanFilter::h(const VectorXd &x) {
    // extract position and velocity
    float px = x(0);
    float py = x(1);
    float vx = x(2);
    float vy = x(3);

    if(fabs(px) < ZERO) {
        std::cout << "Warning: px too small, hack it to : " << ZERO << std::endl;
        px = ZERO;
    }

    float rho = sqrt(px*px + py*py);
    float theta = atan2(py, px);  // arc tangent(2) of y/x, in the interval [-pi,+pi] radians. 
    float rho_dot = (px*vx + py*vy) / rho;

    VectorXd hx = VectorXd(3);
    hx << rho, theta, rho_dot;

    return hx;
}

