%
% run cw1 before run this!!!


%% Compare position between GNSS KF, DR, and corrected after integration
close all;
GNSS_KF_lat = States_NED(:,1)-correct_Deg(1,2); % lat
GNSS_KF_long = States_NED(:,2)-correct_Deg(1,3); % long

DR_lat = L_K-correct_Deg(1,2);
DR_long = Lam_K-correct_Deg(1,3);

% plot corrected position
correct_lat = correct_Deg(:,2)-correct_Deg(1,2); % lat
correct_long = correct_Deg(:,3)-correct_Deg(1,3); % long

hold on
plot(GNSS_KF_lat,GNSS_KF_long, 'b');
plot(DR_lat,DR_long, 'g');
plot(correct_lat,correct_long,'r');
plot(0, 0, 'o')
hold off

legend("GNSS KF", "DR", "Corrected from Integration");
savefig("pos_comparison");
%% Compare velocity North between GNSS KF, DR, and corrected after integration
close all;
t = [1:size(correct_Deg,1)];
v_GNSS = States_NED(:,4);
v_DR = VN;
v_correct = correct_Deg(:,4);

hold on
plot(t,v_GNSS, 'b');
plot(t,v_DR, 'g');
plot(t,v_correct,'r');
hold off

legend("GNSS KF", "DR", "Corrected from Integration");
savefig("vel_north_comparison");
%% Compare velocity East between GNSS KF, DR, and corrected after integration
close all;
t = [1:size(correct_Deg,1)];
v_GNSS = States_NED(:,5);
v_DR = VE;
v_correct = correct_Deg(:,5);

hold on
plot(t,v_GNSS, 'b');
plot(t,v_DR, 'g');
plot(t,v_correct,'r');
hold off

legend("GNSS KF", "DR", "Corrected from Integration");
savefig("vel_east_comparison");
%% plot after multi-epoch positioning
close all;
x = Solutions_NED(:,1); % lat
y = Solutions_NED(:,2); % long
plot(x,y)
savefig("multi_epoch_pos");

%% plot after Kalman Filter with all measurement 
close all;
x = States_out_NED(:,1); % lat
y = States_out_NED(:,2); % long
plot(x,y)
savefig("Kalman_Filter_pos");


%% plot after Kalman Filter GNSS and outlier measurement removed
close all;
x = States_NED(:,1); % lat
y = States_NED(:,2); % long
plot(x,y)
savefig("Kalman_Filter_pos");

%% plot heading compare with original
close all;
t = [1:size(correct_Deg,1)];
plot(t,Heading)
original_heading = D{:,7};

hold on
plot(t,original_heading) 
hold off

legend("Smoothed", "Original");

%% final position | after DR GNSS integration
close all;
num = [1:size(correct_Deg,1)];
% plot corrected position
x = correct_Deg(:,2); % lat
y = correct_Deg(:,3); % long

plot(x,y)
savefig("corrected_pos");

%% final velocity | after DR GNSS integration
close all;
num = [1:size(correct_Deg,1)];
% plot corrected velocity
v_north = correct_Deg(:,4);
v_east = correct_Deg(:,5);

plot(num,v_north);
savefig("corrected_vel_north");
plot(num,v_east);
savefig("corrected_vel_east");


