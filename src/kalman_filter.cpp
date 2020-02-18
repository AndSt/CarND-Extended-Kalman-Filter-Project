#include "kalman_filter.h"
#include "tools.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
  //std::cout << "predict" << std::endl;
   x_ = F_ * x_;
   MatrixXd Ft = F_.transpose();
   P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

  VectorXd y = z - H_ * x_;
  CommonUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  //std::cout << "updateEKF" << std::endl;

  float px= x_[0];
  float py= x_[1];
  float px_dot= x_[2];
  float py_dot= x_[3];

  float rho = sqrt(px*px + py*py);
  if (fabs(rho) < 0.0001){
    rho = 0.0001;
  }
  float theta = atan2(py,px);
  float rho_dot = (px*px_dot + py*py_dot) / rho ;
  VectorXd hx(3);
  hx << rho, theta, rho_dot;

  VectorXd y = z - hx;
  y(1) = norm_pi(y(1));
  CommonUpdate(y);
}

void KalmanFilter::CommonUpdate(const VectorXd &y){
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

double KalmanFilter::norm_pi(double num)
{
  while (num > M_PI)
  {
    num -= 2 * M_PI;
  }
  while (num < -M_PI)
  {
    num += 2 * M_PI;
  }
  return num;
}
