#include "FusionEKF.h"

#include <cmath>
#include <iostream>

#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;
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
    previous_timestamp_ = measurement_pack.timestamp_;

    // first measurement
    cout << "Initializing ... " << endl;

    // init x
    VectorXd x(4);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      x << sin(measurement_pack.raw_measurements_[0]), cos(measurement_pack.raw_measurements_[1]), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // init measurement covariance matrix - laser
    MatrixXd R_laser(2, 2);
    R_laser << 0.0225, 0,
               0, 0.0225;
    // init measurement covariance matrix - radar
    MatrixXd R_radar(3, 3);
    R_radar << 0.09, 0, 0,
               0, 0.0009, 0,
               0, 0, 0.09;

    // set ax, ay to 9
    ekf_.Init(x, R_laser, R_radar, 9, 9);

    // done initializing, no need to predict or update
    is_initialized_ = true;

    cout << "Initialized. " << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  double dt_s = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  ekf_.Predict(dt_s);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  previous_timestamp_ = measurement_pack.timestamp_;

  // print the output
  cout << "x_ = " << ekf_.GetCurrentState() << endl;
  cout << "P_ = " << ekf_.GetCurrentCovariance() << endl;
}
