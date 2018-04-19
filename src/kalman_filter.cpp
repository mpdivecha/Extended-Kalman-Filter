#include "kalman_filter.h"
//#include <cmath>
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define PI 3.14159265

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict()
{
    /**
    TODO:
     * predict the state
    */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

inline VectorXd h_func(VectorXd &x_state)
{
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    VectorXd h(3);
    float denom = sqrt(px * px + py * py);
    h << denom, atan2(py, px),
        (px * vx + py * vy) / denom;

    return h;
}

void KalmanFilter::Update(const VectorXd &z)
{
    /**
    TODO:
     * update the state by using Kalman Filter equations
    */
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    // new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
    /**
    TODO:
     * update the state by using Extended Kalman Filter equations
    */
    VectorXd z_pred = h_func(x_);
    VectorXd y = z - z_pred;
    if (y(1) <= -PI)
    {
        y(1) = y(1) + 2 * PI;
    }
    else if (y(1) >= PI)
    {
        y(1) = y(1) - 2 * PI;
    }

    MatrixXd Hj_ = H_;

    MatrixXd Hjt = Hj_.transpose();
    MatrixXd S = Hj_ * P_ * Hjt + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Hjt;
    MatrixXd K = PHt * Si;

    // new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * Hj_) * P_;
}