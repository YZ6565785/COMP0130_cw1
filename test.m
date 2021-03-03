ts = 0.5; % propagation interval (epoch every 0.5s)
% compute the transition matrix
Phi_k_minus_1 = [
    eye(3) ts*eye(3) zeros(3,1) zeros(3,1);
    zeros(3) eye(3) zeros(3,1) zeros(3,1);
    zeros(1,3) zeros(1,3) 1 ts;
    zeros(1,3) zeros(1,3) 0 1
];

Phi_k_minus_1
size(Phi_k_minus_1)

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

Q_k_minus_1
size(Q_k_minus_1)
