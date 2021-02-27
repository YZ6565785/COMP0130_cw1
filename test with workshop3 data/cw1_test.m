Speed_heading = readtable("Workshop3_Speed_Heading.csv");
Define_Constants;

D = array2table(zeros(size(Speed_heading, 1), 7));
D{:,1} = Speed_heading{:,1};
D{:,2} = Speed_heading{:,2};
D{:,3} = Speed_heading{:,2};
D{:,4} = Speed_heading{:,2};
D{:,5} = Speed_heading{:,2};
% D{:,6} = Speed_heading{:,3};
D{:,7} = Speed_heading{:,3};

States_NED = readtable("Workshop3_GNSS_Pos_Vel_NED.csv");
States_NED = States_NED{:,2:end};
[L_K, Lam_K, damped_vn, damped_ve, VN, VE] = cw1_dead_reckoning(States_NED, D);


% [correct_DR] = cw1_GNSS_DR_integration(States_NED, D, L_K, Lam_K, damped_vn, damped_ve, VN, VE);
[correct_DR] = cw1_DR_GNSS_integration(States_NED, D, L_K, Lam_K, damped_vn, damped_ve, VN, VE);
% convert to degree
correct_Deg = [];
for i=1:size(correct_DR,2)
    lat_deg = correct_DR(2,i) * rad_to_deg;
    long_deg = correct_DR(3,i) * rad_to_deg;
    correct_Deg = [
        correct_Deg; 
        correct_DR(1,i) lat_deg long_deg correct_DR(4,i) correct_DR(5,i)];
end
% correct_Deg