%
% run cw1 before run this!!!


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


%% plot after Kalman Filter and outlier measurement removed
close all;
x = States_NED(:,1); % lat
y = States_NED(:,2); % long
plot(x,y)
savefig("Kalman_Filter_pos");

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

