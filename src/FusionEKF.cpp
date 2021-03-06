#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
* Constructor
*/
FusionEKF::FusionEKF()
{
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initialize matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    // measurment covariance matrix - laser
    R_laser_ << 0.0225, 0,
        0, 0.0225;

    // measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

    /**
    TODO:
     * Finish initializing the FusionEKF.
     * Set the process and measurement noises
    */
    //state covariance matrix P
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;

    H_laser_ << 1, 0, 0, 0,
        0, 1, 0, 0;
}

/**
* Destructor
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
    /*****************************************************************************
     * Initialization
     *****************************************************************************/
    if (!is_initialized_)
    {
        /**
         TODO:
          * Intialize the state ekf_.x_ with the first measurement
          * Create the covariance matrix
          * Remember: you'll need to convert radar from polar to cartesian coordinates
        */
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        previous_timestamp_ = measurement_pack.timestamp_;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
        {
            /**
            Convert radar from polar to cartesian coordinates and initialize state
            */
            float rho = measurement_pack.raw_measurements_[0];
            float phi = measurement_pack.raw_measurements_[1];
            //float rho_dot = measurement_pack(2);

            float px = rho * cos(phi);
            float py = rho * sin(phi);

            ekf_.x_ = VectorXd(4);
            ekf_.x_ << px, py, 0, 0;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
        {
            /**
            Initialize state
            */
            ekf_.x_ = VectorXd(4);
            //set the state with the initial location and zero velocity
            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
        }

        cout << "Initialized" << endl;
        // done initializing, no need to predit or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
        * Prediction
        *****************************************************************************/

    /**
        TODO:
         * Update the state transition matrix F according to the new elapsed time.
          - Time is measured in seconds
         * Update the process noise covariance matrix
         * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
        */

    //compute the time elapsed between the current and previous measurements
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;

    // Update the state transition matrix F
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, dt, 0,
        0, 1, 0, dt,
        0, 0, 1, 0,
        0, 0, 0, 1;

    // Update the process covariance matrix Q
    float dt2 = dt * dt;
    float dt3 = dt2 * dt;
    float dt4 = dt3 * dt;
    float noise_ax = 9;
    float noise_ay = 9;

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt4 / 4 * noise_ax, 0, dt3 / 2 * noise_ax, 0,
        0, dt4 / 4 * noise_ay, 0, dt3 / 2 * noise_ay,
        dt3 / 2 * noise_ax, 0, dt2 * noise_ax, 0,
        0, dt3 / 2 * noise_ay, 0, dt2 * noise_ay;

    ekf_.Predict();

    /*****************************************************************************
         * Update
         ****************************************************************************/

    /**
         TODO:
          * Use the sensor type to perform the update step
          * Update the state and covariance matrices
        */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
        // Radar updates
        // Radar update needs the Jacobian of the H matrix
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.H_ = Hj_;
        // Measurement covariance matrix for Radar
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }
    else
    {
        // Laser updates
        ekf_.H_ = H_laser_;
        // Measurement covariance matrix for Laser
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}