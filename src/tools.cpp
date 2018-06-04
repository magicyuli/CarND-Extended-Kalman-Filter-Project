#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  // ... your code here
  if (estimations.empty() || estimations.size() != ground_truth.size()) {
     throw std::invalid_argument("Invalid dimensions of inputs");
  }

  VectorXd sumSq = VectorXd::Zero(estimations[0].size());
  for(int i = 0; i < estimations.size(); i++) {
    VectorXd diff = estimations[i] - ground_truth[i];
    sumSq = sumSq.array() + diff.array() * diff.array();
  }
  return (sumSq / estimations.size()).array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //check division by zero
  if (px == 0 && py == 0) {
      std::cout << "Division by 0!" << std::endl;
      return MatrixXd::Zero(3, 4);
  }
  
  MatrixXd J(3,4);
  //compute the Jacobian matrix
  float sumSq = px * px + py * py;
  float sqrtSumSq = sqrt(sumSq);
  float sqrtSumSq_3 = sumSq * sqrtSumSq;
  
  J << px/sqrtSumSq, py/sqrtSumSq, 0, 0,
        -py/sumSq, px/sumSq, 0, 0,
        py*(py*vx - px*vy)/sqrtSumSq_3, px*(px*vy - py*vx)/sqrtSumSq_3, px/sqrtSumSq, py/sqrtSumSq;
  
  return J;
}
