#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
    std_a_ = 0.3;
    
    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.3;
    
    //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
    //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
    
    /**
     TODO:
     
     Complete the initialization. See ukf.h for other member properties.
     
     Hint: one or more values initialized above might be wildly off...
     */
    is_initialized_ = false;
    n_x_ = 5;
    n_aug_ = 7;
    lambda_ = 3 - n_aug_;
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    weights_ = VectorXd(2 * n_aug_ + 1);
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1;i++){
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }
    
    //vector<VectorXd> NIS_lidar;
    //vector<VectorXd> NIS_radar;
    
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
    if (!is_initialized_){
        cout << "UKF: " << endl;
        x_ <<   1,1,1,1,1;
        
        P_ <<   1,0,0,0,0,
                0,1,0,0,0,
                0,0,1,0,0,
                0,0,0,1,0,
                0,0,0,0,1;
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            double px = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
            double py = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
            
            x_ << px, py, 0, 0, 0;
            
            time_us_ = meas_package.timestamp_;
        }
        
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0 ,0 ,0;
            time_us_ = meas_package.timestamp_;
        }
        
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
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
    
    // 1) starting with generating augmented sigma points from the modified x state and P covariance
    VectorXd x_aug = VectorXd(n_aug_);
    MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
    
    x_aug.head(n_x_) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;
    
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_,n_x_) = P_;
    MatrixXd noise_matrix = MatrixXd(n_aug_ - n_x_, n_aug_ - n_x_);
    noise_matrix << std_a_*std_a_, 0,
                    0, std_yawdd_*std_yawdd_;
    P_aug.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) = noise_matrix;
    MatrixXd L = P_aug.llt().matrixL(); // square root of augmented P
    
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }
    
    // 2) Predicting sigma points after delta_t
    for(int i=0; i<(2 * n_aug_ + 1); i++){
        VectorXd current_sig = Xsig_aug.col(i);
        VectorXd x_k = VectorXd(n_x_);
        x_k = current_sig.head(n_x_);
        VectorXd noise_part = VectorXd(n_x_);
        noise_part <<   0.5*pow(delta_t,2)*cos(current_sig[3])*current_sig[5],
                        0.5*pow(delta_t,2)*sin(current_sig[3])*current_sig[5],
                        delta_t*current_sig[5],
                        0.5*pow(delta_t,2)*current_sig[6],
                        delta_t*current_sig[6];
        
        VectorXd integral_part = VectorXd(n_x_);
        if(current_sig[4] <= abs(0.001)){
            integral_part <<    current_sig[2]*cos(current_sig[3])*delta_t,
                                current_sig[2]*sin(current_sig[3])*delta_t,
                                0,
                                current_sig[4]*delta_t,
                                0;
        }
        else{
            integral_part <<    current_sig[2] / current_sig[4]*(sin(current_sig[3]+current_sig[4]*delta_t)-sin(current_sig[3])),
                                current_sig[2] / current_sig[4]*(-cos(current_sig[3]+current_sig[4]*delta_t)+cos(current_sig[3])),
                                0,
                                current_sig[4]*delta_t,
                                0;
        }
        Xsig_pred_.col(i) = x_k + integral_part + noise_part;
    }
    
    // 3) Predicting x state and P covariance using weights
    for (int i = 0; i<2*n_aug_+1; i++){
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }
    
    for(int i = 0; i<2*n_aug_+1; i++){
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
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
    
    // 1) Find sigma point in the measurement space Zsig(2*15) by transforming Xsig_pred_
    int n_z = 2;
    MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
    for (int i=0 ; i<2*n_aug_+1; i++){
        Zsig(0,i) = Xsig_pred_(0,i);
        Zsig(1,i) = Xsig_pred_(1,i);
    }
    
    // 2) Find measurement mean z(2*1) and convariance S(2*2) using Zsig with adding noise R(2*2)
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0 ; i<2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i=0 ; i<2*n_aug_+1; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    //add measurement noise covariance matrix
    MatrixXd R_laser_ = MatrixXd(2,2);
    R_laser_ <<     std_laspx_*std_laspx_, 0,
                    0, std_laspy_*std_laspy_;
    
    S = S + R_laser_;
    
    // 3) update x state and P covariance
    //create matrix for cross correlation Tc(5*2)
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i=0 ; i<2*n_aug_+1; i++) {
        // measurement difference
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    VectorXd z = VectorXd(n_z); // incoming Lidar measurement
    z <<    meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1];
    
    //residual
    VectorXd z_diff = z - z_pred;
    
    //update state mean and covariance matrix
    x_= x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    //double epsilon = z_diff.transpose() * S.inverse() * z_diff;
    //NIS_lidar.push_back(epsilon);
    double NIS_lidar = z_diff.transpose() * S.inverse() * z_diff;
    //NIS_lidar.push_back(epsilon);
    
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
    
    // 1) Find sigma point in the measurement space Zsig(3*15) by transforming Xsig_pred_
    int n_z = 3;
    MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
    for(int i=0 ; i<2*n_aug_+1; i++){
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
        
        //angle normalization
        while (Zsig(1,i)> M_PI) Zsig(1,i)-=2.*M_PI;
        while (Zsig(1,i)<-M_PI) Zsig(1,i)+=2.*M_PI;
    }
    
    // 2) Find measurement mean z(3*1) and convariance S(3*3) using Zsig with adding noise R(3*3)
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0 ; i<2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i=0 ; i<2*n_aug_+1; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    //add measurement noise covariance matrix
    MatrixXd R_radar_ = MatrixXd(3,3);
    R_radar_ <<     std_radr_*std_radr_, 0, 0,
                    0, std_radphi_*std_radphi_, 0,
                    0, 0,std_radrd_*std_radrd_;
    S = S + R_radar_;
    
    // 3) update x state and P covariance
    //create matrix for cross correlation Tc(5*3)
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i=0 ; i<2*n_aug_+1; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    VectorXd z = VectorXd(n_z); // incoming radar measurement
    z <<    meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            meas_package.raw_measurements_[2];
    
    //residual
    VectorXd z_diff = z - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    //update state mean and covariance matrix
    x_= x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    //double epsilon = z_diff.transpose() * S.inverse() * z_diff;
    //NIS_radar.push_back(epsilon);
    
    double NIS_radar = z_diff.transpose() * S.inverse() * z_diff;
    //NIS_radar.push_back(epsilon);
}
