#include "kalman_filter.h"

#include <cmath>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(const VectorXd &x_in, const MatrixXd &R_laser_in, const MatrixXd &R_radar_in,
                        const double &ax_in, const double &ay_in) {
  x_ = x_in;
  P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;
  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;
  H_ = MatrixXd(2, 4);
  H_ << 1, 0, 0, 0,
        0, 1, 0, 0;
  R_laser_ = R_laser_in;
  R_radar_ = R_radar_in;
  ax_ = ax_in;
  ay_ = ay_in;
}

void KalmanFilter::Predict(const double &dt_s) {
  UpdateFAndQ(dt_s);
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  InternalUpdate(y, H_, R_laser_);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];
  double rho = sqrt(px * px + py * py);

  VectorXd hx(3);
  hx << rho, atan2(py, px), (px * vx + py * vy) / rho;

  MatrixXd Hj = tools_.CalculateJacobian(x_);
  
  VectorXd my_z = z;
  AdjustRho(my_z, hx);
  InternalUpdate(my_z - hx, Hj, R_radar_);
}

void KalmanFilter::InternalUpdate(const VectorXd &y, const MatrixXd &H, const MatrixXd &R) {
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd K = P_ * Ht * S.inverse();

  //new estimate
  x_ += (K * y);

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;
}

void KalmanFilter::UpdateFAndQ(const double &dt_s) {
  float dt_2 = dt_s * dt_s;
  float dt_3 = dt_2 * dt_s;
  float dt_4 = dt_3 * dt_s;

  //Modify the F matrix so that the time is integrated
  F_(0, 2) = dt_s;
  F_(1, 3) = dt_s;
  
  Q_ = MatrixXd(4, 4);
  Q_ << dt_4/4*ax_, 0, dt_3/2*ax_, 0,
        0, dt_4/4*ay_, 0, dt_3/2*ay_,
        dt_3/2*ax_, 0, dt_2*ax_, 0,
        0, dt_3/2*ay_, 0, dt_2*ay_;
}

static const double PI  = 3.141592653589793238463;
static const double TWO_PI  = 2 * 3.141592653589793238463;

inline void KalmanFilter::AdjustRho(VectorXd &z, const VectorXd &hx) {
  // adjust rho to [-pi, pi]
  while (z(1) > PI) {
    z(1) -= TWO_PI;
  }
  while (z(1) < -PI) {
    z(1) += TWO_PI;
  }

  // adjust rho so that it reflects the real difference from the estimated rho.
  // example of not reflecting the real difference: 3.13 and -3.13, where the real
  // difference is 0.02, instead of 6.26.
  // only need to do this step if the signs of the two rho's are different.
  if (z(1) * hx(1) < 0) {
    double tmp = abs(z(1) - hx(1)) < abs(z(1) + TWO_PI - hx(1)) ? z(1) : z(1) + TWO_PI;
    z(1) = abs(tmp - hx(1)) < abs(z(1) - TWO_PI - hx(1)) ?tmp : z(1) - TWO_PI;
  }
}
