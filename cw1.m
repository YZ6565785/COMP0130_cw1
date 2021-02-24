

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
% [x_plus] = cw1_single_epoch_positioning(Ranges);
% 
% % print the output
% [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_plus(1:3),[0;0;0]);
% lat_deg = L_b * rad_to_deg;
% long_deg = lambda_b * rad_to_deg;
% disp("Answer:");
% disp("lat: " + num2str(lat_deg,'%.6f'));
% disp("long: " + num2str(long_deg,'%.6f'));
% disp("height: " + h_b);

% ================================================= Multi-epoch Positioning
% r_ea_e = [L_b; lambda_b; h_b];
% [Positions, outliers] = cw1_multi_epoch_positioning(Ranges, r_ea_e);
% Pos_NED = [];
% for i=1:size(Positions,2)
%     [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(Positions(1:3,i),[0;0;0]);
%     lat_deg = L_b * rad_to_deg;
%     long_deg = lambda_b * rad_to_deg;
%     Pos_NED = [Pos_NED; lat_deg long_deg h_b];
% end
% Pos_NED
% 
% % print the outlier detected
% for i=1:size(outliers, 1)
%     disp("Found outlier at time: ["+outliers(i,1)+"], iter: ["+...
%         outliers(i,2)+"], satellite: ["+outliers(i,3)+"], residual: ["+...
%         outliers(i,4)+"]");
% end


% ================================================= Multi-epoch Positioning
[States] = cw1_GNSS_kalman_filter_multi(Ranges, Rates);




