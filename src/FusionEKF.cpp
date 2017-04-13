#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0ULL;

  // initializing matrices
  F_ = MatrixXd(4, 4);
  Q_ = MatrixXd(4, 4);
 
  H_laser_ = MatrixXd(2, 4);
  Hj_      = MatrixXd(3, 4);
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);

  //laser measurement matrix
  H_laser_ << 1, 0, 0, 0,
           0, 1, 0, 0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  //set the acceleration noise components
  noise_ax = 9.f;
  noise_ay = 9.f;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

    /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
      cout << "EKF: " << endl;
      ekf_.x_ = VectorXd(4);
      ekf_.P_ = MatrixXd(4, 4);
      //ekf_.x_ << 1, 1, 1, 1;
      //state covariance matrix P
      ekf_.P_ << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1000, 0,
         0, 0, 0, 1000;

      if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        /**
          Convert radar from polar to cartesian coordinates and initialize state.
         */
        double rho = measurement_pack.raw_measurements_[0]; // Range - radial distance from origin
        double phi = measurement_pack.raw_measurements_[1]; // Bearing - angle between rho and x
        double rho_dot = measurement_pack.raw_measurements_[2]; // Radial Velocity - change of p (range rate)
        double x = rho * cos(phi);
        double y = rho * sin(phi);
        double vx = rho_dot * cos(phi);
        double vy = rho_dot * sin(phi);

        ekf_.x_ << x, y, vx, vy;  // division by zero is handled in filter logic
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        /**
          Initialize state.
         */
        // Ladar has no velocity data
        ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    
    previous_timestamp_ = measurement_pack.timestamp_;
    
    ekf_.Init(ekf_.x_, ekf_.P_, F_, H_laser_, R_laser_, Q_);  // division by zero is handled in filter logic
    // done initializing, no need to predict or update
    is_initialized_ = true;
    cout << "init x_ = " << ekf_.x_ << endl;
    cout << "init P_ = " << ekf_.P_ << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  //compute the time elapsed between the current and previous measurements
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  //the initial transition matrix F_
  F_ << 1, 0, dt, 0,
     0, 1, 0, dt,
     0, 0, 1, 0,
     0, 0, 0, 1;

  // Update the process noise covariance matrix.
  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;
  Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
     0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
     dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
     0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  // update F, Q
  ekf_.Init(ekf_.x_, ekf_.P_, F_, ekf_.H_, ekf_.R_, Q_);
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar updates
      VectorXd z = VectorXd(3);
      z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1],
        measurement_pack.raw_measurements_[2];

      Hj_ = Tools::CalculateJacobian(ekf_.x_);
      //std::cout << "Hj_ " << Hj_ << std::endl;
      // update H, R
      ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);
      ekf_.UpdateEKF(z);

  } else {
      // Laser updates
      VectorXd z = VectorXd(2);
      z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
      // update H, R
      ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, H_laser_, R_laser_, ekf_.Q_);
      ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
