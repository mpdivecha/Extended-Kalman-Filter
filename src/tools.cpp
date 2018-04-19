#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    /**
    TODO:
     * Calculate the RMSE here
    */
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    // check the validity of the following inputs:
    // * the estimation vector size should not be zero
    // * the estimation vector size should equal ground truth vector size
    if (estimations.size() != ground_truth.size()
            || estimations.size() == 0) {
                cout << "Invalid estimation or ground_truth data" << endl;
                return rmse;
    }

    // accumulate sqaured residuals
    for (unsigned int i = 0; i < estimations.size(); ++i) {

        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficeint-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }

    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the sqaured root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    /**
    TODO:
     * Calculate a Jacobian here
    */
    MatrixXd Hj(3, 4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    // check for division by zero
    if (px == 0 && py == 0) {
        cout << "Division by zero" << endl;
        return Hj;
    }

    //cout << "Computing needed variables" << endl;
    // compute the Jacobian matrix
    float denom = px*px + py*py;
    float denomSqrt = sqrt(denom);
    float denom3Sqrt = sqrt(denom*denom*denom);

    Hj << px/denomSqrt, py/denomSqrt, 0, 0,
          -py/denom   , px/denom    , 0, 0,
          py*(vx*py - vy*px)/denom3Sqrt, px*(vy*px - vx*py)/denom3Sqrt, px/denomSqrt, py/denomSqrt;
    
    return Hj;
}