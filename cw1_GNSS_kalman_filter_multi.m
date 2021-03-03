function [States] = cw1_GNSS_kalman_filter_multi(Ranges, Rates, state_init, P_matrix)
% cw1_GNSS_kalman_filter_multi() - compute the position at multiple epochs 
% using the all rows of the pseudo-range measurements. Starting by an 
% initial position at the centre of the Earth.
%
% Reference:
% COMP0130: ROBOT VISION AND NAVIGATION
% Workshop 1: Task 2
%
% This function created 22/02/2021 by Yuhang Zhang
%
% Inputs:
%   Ranges           Pesudo-ranges table read from files.
%   state_init       initial state (in ECEF).
%                    row 1: position x
%                    row 2: position y
%                    row 3: position z
%                    row 4: velocity x
%                    row 5: velocity y
%                    row 6: velocity z
%                    row 7: clock offset
%                    row 8: clock drift
%   P_matrix         initial error covariance matrix.
%
% Outputs:
%   States           Computed positions and velocities for every epochs 
%                    with outliers removed. (ECEF).

% Copyright 2021, Yuhang Zhang
% License: BSD; see license.txt for details

% =========================================================================
% Begins
States = [];

% load constants
Define_Constants;

% time list
times = Ranges{2:end,1};
% satellite number list
numbers = Ranges{1, 2:end};

% initialize the Kalman filter state vetor estimate 
% & error covariance matrix
% [x_est,P_matrix] = Initialise_GNSS_KF;
x_0_plus = state_init;
P_0_plus = P_matrix;


ts = 0.5; % propagation interval (epoch every 0.5s)
% compute the transition matrix
Phi_k_minus_1 = [
    eye(3) ts*eye(3) zeros(3,1) zeros(3,1);
    zeros(3) eye(3) zeros(3,1) zeros(3,1);
    zeros(1,3) zeros(1,3) 1 ts;
    zeros(1,3) zeros(1,3) 0 1
];


% power spectral desity (PSD)
S_a_e = 5;
% clock phase PSD
S_c_phi_a = 0.01;
% clock frequency PSD
S_c_f_a = 0.04;

% compute noise covariance matrix
Q_k_minus_1 = [
    S_a_e*ts^3*eye(3)/3 S_a_e*ts^2*eye(3)/2 zeros(3,1) zeros(3,1);
    S_a_e*ts^2*eye(3)/2 S_a_e*ts*eye(3) zeros(3,1) zeros(3,1);
    zeros(1,3) zeros(1,3) S_c_phi_a*ts+S_c_f_a*ts^3/3 S_c_f_a*ts^2/2;
    zeros(1,3) zeros(1,3) S_c_f_a*ts^2/2 S_c_f_a*ts;
];


for epoch=1:size(times,1)
% for epoch=1:2
    % propagate the state estimates:
    x_k_minus = Phi_k_minus_1*x_0_plus;
%     x_k_minus
    
    % propagate the error covariance matrix
    P_k_minus = Phi_k_minus_1*P_0_plus*Phi_k_minus_1'+Q_k_minus_1;
    %     P_k_minus

    % initialise postion and velocity of each satellite
    r_ej_e = zeros(3, size(numbers, 2));
    v_ej_e = zeros(3, size(numbers, 2));
    
    % initialise the range rates
    r_aj_minus = zeros(size(numbers,2),1);

    % initialise the line-of-sight unit vector
    u_aj_e = zeros(3, size(numbers,2));
    
    % initialise the range rates
    r_dot_aj_minus = zeros(size(numbers,2), 1);
    
    for i=1:size(numbers,2)
        % predict the ranges from approximated user position
        [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(times(epoch), numbers(i));
        r_ej_e(:, i) = sat_r_es_e;
        v_ej_e(:, i) = sat_v_es_e;
        
        f = eye(3) * r_ej_e(:,i) - x_k_minus(1:3);
        r_a = sqrt(f'*f);
        
        for itr=1:5
            C_e_I = [
                1 omega_ie*r_a/c 0;
                -omega_ie*r_a/c 1 0;
                0 0 1
            ];

            f = C_e_I * r_ej_e(:,i) - x_k_minus(1:3);
            r_a = sqrt(f'*f);
        end
        r_aj_minus(i) = r_a;
        C_e_I = [
            1 omega_ie*r_aj_minus(i)/c 0;
            -omega_ie*r_aj_minus(i)/c 1 0;
            0 0 1
        ];
        % compute the line-of-sight unit vector
        u_aj_e(:,i) = (C_e_I*r_ej_e(:,i) - x_k_minus(1:3)) / r_aj_minus(i);
        
        % predict the range rates from approximated user position
        r_dot_aj_minus(i) = u_aj_e(:,i)' * (...
            C_e_I*(v_ej_e(:,i)+Omega_ie*r_ej_e(:,i))- ...
                (x_k_minus(4:6)+Omega_ie*x_k_minus(1:3)) ...
        );
    end
    
    % compute the measurement matrix 
    H_k = [
        -u_aj_e' zeros(size(numbers,2),3) ...
        ones(size(numbers,2),1) zeros(size(numbers,2),1);
        zeros(size(numbers,2),3) -u_aj_e' ...
        zeros(size(numbers,2),1) ones(size(numbers,2),1)
    ];

    % compute the measurement noise covariance matrix
    sigma_p = 10;
    sigma_r = 0.05;
    R_k = eye(size(numbers,2)*2);
    for i=1:size(R_k,1)
        if i<=size(numbers,2)
            R_k(i,i) = sigma_p^2;
        else
            R_k(i,i) = sigma_r^2;
        end
    end
%     R_k
    % compute the Kalman gain matrix
    K_k = P_k_minus*H_k'/(H_k*P_k_minus*H_k'+R_k);

    % compute the measurement innovation vector 
    ranges = Ranges{epoch+1, 2:end}';
    rates = Rates{epoch+1, 2:end}';
    delta_z = [
        ranges-r_aj_minus-x_k_minus(7);
        rates-r_dot_aj_minus-x_k_minus(8)
    ];
%     delta_z

    % update the state estimates
    x_k_plus = x_k_minus + K_k*delta_z;
    %     x_k_plus
    
    % update the error covariance matrix
    P_k_plus = (eye(8)-K_k*H_k)*P_k_minus;
%     P_k_plus
    
    % store the state at each epoch
%     [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_k_plus(1:3),x_k_plus(4:6));
%     disp("Answer:");
%     disp("latitude: " + L_b * rad_to_deg );
%     disp("longitude: " + lambda_b * rad_to_deg );
%     disp("height: " + h_b );
%     disp(v_eb_n);
%     States = [States; times(epoch) L_b * rad_to_deg lambda_b * rad_to_deg h_b v_eb_n'];
    States = [States x_k_plus(1:6)];
    
    
    % prepare for next epoch
    x_0_plus = x_k_plus;
    P_0_plus = P_k_plus;
end



% Ends