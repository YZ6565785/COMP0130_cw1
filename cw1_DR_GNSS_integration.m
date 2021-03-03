function [Corrected] = cw1_DR_GNSS_integration(States_NED, D, L_K, Lam_K, VN, VE)
% load constants
Define_Constants;

Corrected = [];



% define constants ====================================
% velocity uncertainty
sigma_v = 0.05;
% position uncertainty
sigma_r = 10;
% propagation interval
ts = 0.5;
% GNSS position error standard deviation
sig_Gr = 10;
% GNSS velocity error standard deviation
sig_Gv = 0.05;
% DR velocity error power spectral density (PSD)
S_DR = 0.01;
% =====================================================

% declare arrays
times = D{:,1};

lats_G = States_NED(:,1)*deg_to_rad;
longs_G = States_NED(:,2)*deg_to_rad;
v_n_G = States_NED(:,4);
v_e_G = States_NED(:,5);
heights = States_NED(:,3);

lats_D = L_K*deg_to_rad;
longs_D = Lam_K*deg_to_rad;
v_n_D = VN;
v_e_D = VE;

% init states and matrix
x_0_plus = [0;0;0;0];
lat_0 = lats_G(1);
[R_N,R_E]= Radii_of_curvature(lat_0);
h_0 = heights(1);
P_0_plus = diag([
    sigma_v^2; 
    sigma_v^2; 
    sigma_r^2/(R_N+h_0)^2; 
    sigma_r^2/((R_E+h_0)^2*cos(lat_0)^2)
]);

% initialize for time0
Corrected = [Corrected [
    times(1);
    lats_G(1);
    longs_G(1);
    v_n_G(1);
    v_e_G(1)
]];

for epoch=2:size(times,1)
    h_0 = heights(epoch-1);
%     lat_0 = L(epoch-1)*deg_to_rad;
    
    % transition matrix
    Phi_0 = eye(4);
    Phi_0(3,1) = ts/(R_N+h_0);
    Phi_0(4,2) = ts/((R_E+h_0)*cos(lat_0));
    % Phi_0
    
    % noise covariance matrix
    Q_0 = [
        S_DR*ts 0 1/2*S_DR*ts^2/(R_N+h_0) 0;
        0 S_DR*ts 0 1/2*S_DR*ts^2/((R_E+h_0)*cos(lat_0));
        1/2*S_DR*ts^2/(R_N+h_0) 0 1/3*S_DR*ts^3/(R_N+h_0)^2 0;
        0 1/2*S_DR*ts^2/((R_E+h_0)*cos(lat_0)) 0 1/3*S_DR*ts^3/((R_E+h_0)^2*cos(lat_0)^2);
    ];
    % state estimates
    x_k_minus = Phi_0*x_0_plus;
    %error covariance matrix
    P_k_minus = Phi_0*P_0_plus*Phi_0'+Q_0;
%     P_k_minus

    H_k = zeros(4);
    H_k(3,1) = -1;
    H_k(4,2) = -1;
    H_k(1,3) = -1;
    H_k(2,4) = -1;
    
    % meridian radius and transverse radius
    lat_k = lats_G(epoch);
    [R_N,R_E]= Radii_of_curvature(lat_k);
    h_k = heights(epoch);
    
    % measurement noise covariance matrix
    R_k = diag([
        sig_Gr^2/(R_N+h_k)^2;
        sig_Gr^2/((R_E+h_k)^2*cos(lat_k)^2);
        sig_Gv^2;
        sig_Gv^2;
    ]);
    
    % Kalman gain matrix
    K_k = P_k_minus*H_k'/(H_k*P_k_minus*H_k'+R_k);
    
    % measurement innovation vector
    d_z = [
        (lats_G(epoch)-lats_D(epoch));
        (longs_G(epoch)-longs_D(epoch));
        v_n_G(epoch) - v_n_D(epoch);
        v_e_G(epoch) - v_e_D(epoch);
    ] - H_k*x_k_minus;

    % update state estimates
    x_k_plus = x_k_minus + K_k*d_z;
    % update the error covariance matrix
    P_k_plus = (eye(4)-K_k*H_k)*P_k_minus;
    
    x_0_plus = x_k_plus;
    P_0_plus = P_k_plus;
    
    % compute corrected data
    lat_C = lats_D(epoch) - x_k_plus(3);
    long_C = longs_D(epoch) - x_k_plus(4);
    v_nk_C = v_n_D(epoch) - x_k_plus(1);
    v_ek_C = v_e_D(epoch) - x_k_plus(2);
    
    Corrected = [Corrected [
        times(epoch);
        lat_C;
        long_C;
        v_nk_C;
        v_ek_C
    ]];
    
    % prepare for next epoch
    lat_0 = lat_C;
    
end

