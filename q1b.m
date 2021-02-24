clear All;

disp(" ");
disp("Task 1b ===========");

Define_Constants;

% cartesian ECEF positions of satellites
p = readtable('Workshop1_Pseudo_ranges.csv');
disp("cartesian ECEF positions of satellites:");
r_ej_e_0_list = zeros(8, 3);
for i=2:9
    [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(0, p{1,i});
    r_ej_e_0_list(i-1,:) = sat_r_es_e;
end
disp(r_ej_e_0_list);

lat_deg = -33.8125;
long_deg = 151.1973;
h = 60.6318;

[r_eb_e,v_eb_e] = pv_NED_to_ECEF(lat_deg * deg_to_rad, long_deg * deg_to_rad, h, [0;0;0]);

[x_hat_plus] = iter_prediction(r_ej_e_0_list, p, r_eb_e);

[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_hat_plus(1:3),[0;0;0]);

disp("Answer:");
disp("lat: " + L_b * rad_to_deg);
disp("long: " + lambda_b * rad_to_deg);
disp("height: " + h_b);