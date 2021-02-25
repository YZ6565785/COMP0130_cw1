function [C_DR] = cw1_GNSS_DR_integration(States_NED, D)

% load constants
Define_Constants;

% Dead Reckoning
DR = table2array(D);
% average speed
v_mean = (DR(:,2)+DR(:,3)+DR(:,4)+DR(:,5))/4;
% heading
heading = DR(:,7)*deg_to_rad;
V_n = [];
V_e = [];
% average velocity
for i = 2:length(heading)
    v_n = (1/2)*(cos(heading(i)) + cos(heading(i-1)))*v_mean(i);
    v_e = (1/2)*(sin(heading(i)) + sin(heading(i-1)))*v_mean(i);
    V_n = [V_n; v_n];
    V_e = [V_e; v_e];
end
% time
t = DR(:,1);
% inital geodectic latitude and longitude at time 0
L = States_NED(1,1)*deg_to_rad;
Lam = States_NED(1,2)*deg_to_rad;
% height
h = States_NED(:,3);
% DR latitude and longitude 
L_K = [L];
Lam_K = [Lam];
for i = 2:length(t)
    [R_N, R_E] = Radii_of_curvature(L);
    L = L + V_n(i-1)*(t(i)-t(i-1))/(R_N + h(i));
    Lam = Lam + V_e(i-1)*(t(i)-t(i-1))/((R_E + h(i))*cos(L));
    L_K = [L_K; L];
    Lam_K = [Lam_K; Lam];
end
L_K = L_K*rad_to_deg;
Lam_K = Lam_K*rad_to_deg;

% damped instantaneous DR velocity
damped_vn = v_mean(1)*cos(heading(1));
damped_ve = v_mean(1)*sin(heading(1));
VN = [damped_vn];
VE = [damped_ve];

for i = 1:length(V_n)
    damped_vn = 1.7*V_n(i) - 0.7*damped_vn;
    damped_ve = 1.7*V_e(i) - 0.7*damped_ve;
    VN = [VN; damped_vn];
    VE = [VE; damped_ve];
end

%%Integration
% initial state error
x_plus = [0; 0; 0; 0];
% height, latitude and corresponding
% meridian radius and transverse radius at time 0
h0 = h(1);
L0 = States_NED(1,1)*deg_to_rad;
[R_N0, R_E0] = Radii_of_curvature(L0);
% initial state estimation error covariance matrix
p_plus = [0.05^2 0 0 0;...
    0 0.05^2 0 0;... 
    0 0 100/(R_N0 + h0) 0;...
    0 0 0 100/((R_E0 + h0)^2*(cos(L0)^2))];
% propagation interval
ts = 0.5;
% GNSS position error standard deviation
Gr = 10;
% GNSS velocity error standard deviation
Gv = 0.05;
% DR velocity error power spectral density (PSD)
S = 0.01;
% GNSS-indicated solution
GNSS = States_NED(:, [1,2,4,5]);
GNSS(:,1) = GNSS(:,1)*deg_to_rad;
GNSS(:,2) = GNSS(:,2)*deg_to_rad;
% DR-indicated solution
DR_states = [L_K*deg_to_rad Lam_K*deg_to_rad VN VE];
% latitude
L = DR_states(:,1);
C_DR = [DR_states(1,:)];

for i= 2:length(t)
    % transition matrix
    transition = eye(4);
    % meridian radius and transverse radius
    [R_N, R_E] = Radii_of_curvature(L(i));
    transition(3,1) = ts/(R_N + h(i-1));
    transition(4,2) = ts/((R_E +h(i-1))*cos(L(i-1)));
    % noise covariance matrix
    q_con = [S*ts 0 (1/2)*(S*ts^2/(R_N + h(i-1))) 0;...
        0 S*ts 0 (1/2)*(S*ts^2/((R_N + h(i-1))*cos(L(i-1))));...
        (1/2)*(S*ts^2/(R_N + h(i-1))) 0 (1/3)*(S*ts^3/(R_N + h(i-1))^2) 0;...
        0 (1/2)*(S*ts^2/((R_N + h(i-1))*cos(L(i-1)))) 0 (1/3)*(S*ts^3/(((R_N + h(i-1))^2)*cos(L(i-1))^2))];
    % state estimates
    x_min = transition*x_plus;
    %error covariance matrix
    p_min = transition*p_plus*transpose(transition)+q_con;
    % measurement matrix
    H = [0 0 -1 0; 0 0 0 -1; -1 0 0 0; 0 -1 0 0];
    % measurement noise covariance matrix
    R = [Gr^2/(R_N + h(i))^2 0 0 0;...
        0 Gr^2/(((R_E + h(i))^2)*cos(L(i))^2) 0 0;...
        0 0 Gv^2 0;...
        0 0 0 Gv^2];
    % Kalman gain matrix
    Gain = p_min*transpose(H)*inv(H*p_min*transpose(H)+R);
    % measurement innovation vector
    z_min = transpose(GNSS(i,:) - DR_states(i,:)) - H*x_min;
    % update state estimates
    x_plus = x_min + Gain*z_min;
    % update the error covariance matrix
    p_plus = (eye(4) - Gain*H)*p_min;
    % correct the DR solution
    C = DR_states(i,:) - transpose(x_plus);
    C_DR = [C_DR; C];
end
C_DR(:,1) = C_DR(:,1)*rad_to_deg;
C_DR(:,2) = C_DR(:,2)*rad_to_deg;
end