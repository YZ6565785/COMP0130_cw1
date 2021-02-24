

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
disp("height: " + h_b);
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


% ================================================= Multi-epoch Positioning
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
[States] = cw1_GNSS_kalman_filter_multi(Ranges, Rates, state_init, P_matrix);

% convert to NED
States_NED = [];
for i=1:size(States,2)
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(States(1:3,i),States(4:6,i));
    lat_deg = L_b * rad_to_deg;
    long_deg = lambda_b * rad_to_deg;
    States_NED = [States_NED; lat_deg long_deg h_b v_eb_n'];
end
% States_NED



