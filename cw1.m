

% COMP0130 - Group Z
% Cuinai Yuan, Weizheng Zhang, Yuhang Zhang

clear;
Define_Constants;

% =============================================================== load data
%   1         2         3         4         5           6            7
% time s | speed 1 | speed 2 | speed 3 | speed 4 | gyro angular | heading
D = readtable("Dead_reckoning.csv");

Ranges = readtable("Pseudo_ranges.csv");
Rates = readtable("Pseudo_range_rates.csv");


% ========================= Single-epoch Positioning without Initialisation
[x_plus] = cw1_single_epoch_positioning(Ranges);

% print the output
[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_plus(1:3),[0;0;0]);
lat_deg = L_b * rad_to_deg;
long_deg = lambda_b * rad_to_deg;
disp("Answer:");
disp("lat: " + num2str(lat_deg,'%.6f'));
disp("long: " + num2str(long_deg,'%.6f'));
disp("height: " + num2str(h_b,'%.6f'));
disp("v north: " + v_eb_n(1));
disp("v east: " + v_eb_n(2));
disp("v down: " + v_eb_n(3));

% ================================================= Multi-epoch Positioning
r_ea_e = [L_b; lambda_b; h_b];
[Solutions, outliers] = cw1_multi_epoch_positioning(Ranges, Rates, r_ea_e);
% convert to NED
Solutions_NED = [];
for i=1:size(Solutions,2)
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(Solutions(1:3,i),Solutions(4:6,i));
    lat_deg = L_b * rad_to_deg;
    long_deg = lambda_b * rad_to_deg;
    Solutions_NED = [Solutions_NED; lat_deg long_deg h_b v_eb_n'];
end
% Solutions_NED


% print the outlier detected
for i=1:size(outliers, 1)
    disp("Found outlier at time: ["+outliers(i,1)+"], iter: ["+...
        outliers(i,2)+"], satellite: ["+outliers(i,3)+"], residual: ["+...
        outliers(i,4)+"]");
end


% ======================== GNSS Kalman Filter outlier with all measurements
state_init = Solutions(:,1);
P_matrix =  zeros(8);
P_matrix(1,1) = 100;
P_matrix(2,2) = 100;
P_matrix(3,3) = 100;
P_matrix(4,4) = 0.05^2;
P_matrix(5,5) = 0.05^2;
P_matrix(6,6) = 0.05^2;
P_matrix(7,7) = 100000;
P_matrix(8,8) = 200;

[States_out] = cw1_GNSS_kalman_filter_multi(Ranges, Rates, state_init, P_matrix);

% convert to NED
States_out_NED = [];
for i=1:size(States_out,2)
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(States_out(1:3,i),States_out(4:6,i));
    lat_deg = L_b * rad_to_deg;
    long_deg = lambda_b * rad_to_deg;
    States_out_NED = [States_out_NED; lat_deg long_deg h_b v_eb_n'];
end
% ========================== GNSS Kalman Filter outlier measurement removed
% remove measurement which has largest residuals
ind_remove = outliers(1,5);
Ranges(:,ind_remove+1) = [];
Rates(:,ind_remove+1) = [];

[States] = cw1_GNSS_kalman_filter_multi(Ranges, Rates, state_init, P_matrix);

% convert to NED
States_NED = [];
for i=1:size(States,2)
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(States(1:3,i),States(4:6,i));
    lat_deg = L_b * rad_to_deg;
    long_deg = lambda_b * rad_to_deg;
    States_NED = [States_NED; lat_deg long_deg h_b v_eb_n'];
end


% ===================================================== DR/GNSS Integration
[L_K, Lam_K, VN, VE] = cw1_dead_reckoning(States_NED, D);
% [correct_DR] =  cw1_GNSS_DR_integration(States_NED, D, L_K, Lam_K, damped_vn, damped_ve, VN, VE);
% new integration function defined
[correct_DR] = cw1_DR_GNSS_integration(States_NED, D, L_K, Lam_K, VN, VE);
% convert to degree
correct_Deg = []; 
for i=1:size(correct_DR,2)
    lat_deg = correct_DR(2,i) * rad_to_deg;
    long_deg = correct_DR(3,i) * rad_to_deg;
    correct_Deg = [
        correct_Deg; 
        correct_DR(1,i) lat_deg long_deg correct_DR(4,i) correct_DR(5,i)];
end



