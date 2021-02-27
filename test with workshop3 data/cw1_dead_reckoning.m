function [L_K, Lam_K, damped_vn, damped_ve, VN, VE] = cw1_dead_reckoning(States_NED, D)


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