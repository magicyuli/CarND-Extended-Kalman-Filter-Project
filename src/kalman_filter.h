#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include "Eigen/Dense"

#include "tools.h"

class KalmanFilter {
public:

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
   * @param x_in Initial state
   * @param R_laser_in Measurement covariance matrix
   * @param R_radar_in Measurement covariance matrix
   * @param ax_in Process noise
   * @param ay_in Process noise
   */
  void Init(const Eigen::VectorXd &x_in, const Eigen::MatrixXd &R_laser_in,
            const Eigen::MatrixXd &R_radar_in, const double &ax_in, const double &ay_in);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param dt_s Time between k and k+1 in s
   */
  void Predict(const double &dt_s);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z);

  Eigen::VectorXd GetCurrentState() const {
    return x_;
  }

  Eigen::MatrixXd GetCurrentCovariance() const {
    return P_;
  }

private:

  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state transition matrix
  Eigen::MatrixXd F_;

  // process covariance matrix
  Eigen::MatrixXd Q_;

  // measurement matrix
  Eigen::MatrixXd H_;

  // measurement covariance matrices
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;

  // process noises
  double ax_;
  double ay_;

  Tools tools_;

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param y Result of z - h(x') or z - H*x'
   * @param H H_ or Hj
   * @param R R_laser_ or R_radar_
   */
  void InternalUpdate(const Eigen::VectorXd &y, const Eigen::MatrixXd &H, const Eigen::MatrixXd &R);

  void UpdateFAndQ(const double &dt_s);

  /**
   * Adjust rho of the radar measurement to a reasonable range
   * @param z Measurement
   * @param hx h(x')
   */
  inline void AdjustRho(VectorXd &z, const VectorXd &hx);
};

#endif /* KALMAN_FILTER_H_ */
