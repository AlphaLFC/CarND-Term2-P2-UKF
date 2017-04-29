#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.9;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.55;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  /**
   * ANS by AlphaLFC
   */

  // Initialize first
  is_initialized_ = false;

  // time stamp
  time_us_ = 0;

  // State dim
  n_x_ = 5;

  // Augmented state dim
  n_aug_ = 7;

  // Lambda
  lambda_ = 3 - n_aug_;

  // Initialize weights
  weights_ = VectorXd(2*n_aug_+1);
  double weights_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weights_0;
  double weights_i = 0.5/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {
    weights_(i) = weights_i;
  }

  // Covariance initialization
  x_.fill(0);
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // Predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // Measurement deviation matrices
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_,     0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,   std_radrd_*std_radrd_;

  // NIS
  NIS_laser_ = 0;
  NIS_radar_ = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /**
   * ANS by AlphaLFC
   */

  // Initialization
  if (!is_initialized_) {
    cout << "UKF initializing..." << endl;

    // Laser
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }

    // Radar
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
    }

    if (fabs(x_(0)) < 0.001 and fabs(x_(1)) < 0.001) {
      x_(0) = 0.001;
      x_(1) = 0.001;
    }

    // Initialization done
    is_initialized_ = true;
    cout << "UKF initialized." << endl;

    // Time stamp
    time_us_ = meas_package.timestamp_;
    return;
  }


  // Prediction

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  
  Prediction(delta_t);

  // Update
  if (use_laser_ and (meas_package.sensor_type_ == MeasurementPackage::LASER)) {
    UpdateLidar(meas_package);
  }

  if (use_radar_ and (meas_package.sensor_type_ == MeasurementPackage::RADAR)) {
    UpdateRadar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /**
   * ANS by AlphaLFC
   */

  // Augmented sigma points
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  MatrixXd A_aug = P_aug.llt().matrixL();
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  Xsig_aug.col(0) = x_aug;
  for (int i=0; i<n_aug_; i++) {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * A_aug.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A_aug.col(i);
  }

  // Sigma point prediction
  for (int i=0; i<2*n_aug_+1; i++) {
    double px       = Xsig_aug(0, i);
    double py       = Xsig_aug(1, i);
    double v        = Xsig_aug(2, i);
    double yaw      = Xsig_aug(3, i);
    double yawd     = Xsig_aug(4, i);
    double nu_a     = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    double c1 = sin(yaw);
    double c2 = cos(yaw);
    double c3 = 0.5 * delta_t * delta_t;
    if (fabs(yawd) > 0.001) {
      double c4 = v / yawd;
      double c5 = yawd * delta_t;
      double c6 = yaw + c5;
      px   += c4*( sin(c6)-c1) + c3*c2*nu_a;
      py   += c4*(-cos(c6)+c2) + c3*c1*nu_a;
      v    += delta_t*nu_a;
      yaw  += c5 + c3*nu_yawdd;
      yawd += delta_t*nu_yawdd;
    }
    else {
      double c7 = v * delta_t;
      px   += c2*c7 + c3*c2*nu_a;
      py   += c1*c7 + c3*c1*nu_a;
      v    += delta_t*nu_a;
      yaw  += c3*nu_yawdd;
      yawd += delta_t*nu_yawdd;
    }

    Xsig_pred_(0, i) = px;
    Xsig_pred_(1, i) = py;
    Xsig_pred_(2, i) = v;
    Xsig_pred_(3, i) = yaw;
    Xsig_pred_(4, i) = yawd;
  }

  // Predict mean state and covariance
  x_.fill(0);
  P_.fill(0);

  for (int i=0; i<2*n_aug_+1; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  for (int i=0; i<2*n_aug_+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3) >  M_PI) x_diff(3) -= 2*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2*M_PI;
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  /**
   * ANS by AlphaLFC
   */

  // Measurement dimension
  int n_z = 2;

  // raw measurement
  VectorXd z = meas_package.raw_measurements_;

  // predict measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0);

  // sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  for (int i=0; i<2*n_aug_+1; i++) {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }

  // calculate mean predict meas
  for (int i=0; i<2*n_aug_+1; i++)
    z_pred = z_pred + weights_(i) * Zsig.col(i);

  // measurement covariance
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0);
  S += R_laser_;
  for (int i=0; i<2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // cross correlation
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0);
  for (int i=0; i<2*n_aug_+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3) >  M_PI) x_diff(3) -= 2*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2*M_PI;

    VectorXd z_diff = Zsig.col(i) - z_pred;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain
  MatrixXd K = Tc * S.inverse();

  // update
  VectorXd err = z - z_pred;

  x_ += K * err;
  P_ += -K * S * K.transpose();

  NIS_laser_ = err.transpose() * S.inverse() * err;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  /**
   * ANS by AlphaLFC
   */

  // Measurement dimension
  int n_z = 3;

  // raw measurement
  VectorXd z = meas_package.raw_measurements_;

  // predict measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0);

  // sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  for (int i=0; i<2*n_aug_+1; i++) {
    double px   = Xsig_pred_(0, i);
    double py   = Xsig_pred_(1, i);
    double v    = Xsig_pred_(2, i);
    double yaw  = Xsig_pred_(3, i);
    double yawd = Xsig_pred_(4, i);

    double r = sqrt(px*px + py*py);
    double phi = atan2(py, px);
    double r_dot;

    if (fabs(r) > 0.001)
      r_dot = (px*cos(yaw)*v + py*sin(yaw)*v) / r;
    else
      r_dot = 0;

    Zsig(0, i) = r;
    Zsig(1, i) = phi;
    Zsig(2, i) = r_dot;
  }

  // calculate mean predict meas
  for (int i=0; i<2*n_aug_+1; ++i)
    z_pred = z_pred + weights_(i) * Zsig.col(i);

  // measurement covariance
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0);
  S += R_radar_;
  for (int i=0; i<2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1) >  M_PI) z_diff(1) -= 2*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2*M_PI;
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // cross correlation
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0);
  for (int i=0; i<2*n_aug_+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3) >  M_PI) x_diff(3) -= 2*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2*M_PI;

    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1) >  M_PI) z_diff(1) -= 2*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2*M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain
  MatrixXd K = Tc * S.inverse();

  // update
  VectorXd err = z - z_pred;
  while (err(1) >  M_PI) err(1) -= 2*M_PI;
  while (err(1) < -M_PI) err(1) += 2*M_PI;

  x_ += K * err;
  P_ += -K * S * K.transpose();

  NIS_radar_ = err.transpose() * S.inverse() * err;
}
