%
% run cw1 before run this!!!

close all;
num = [1:size(correct_Deg,1)];

% plot corrected position
x = correct_Deg(:,2);
y = correct_Deg(:,3);
% y = [1:size(Solutions_NED,1)];
plot(x,y)
savefig("corrected_pos");


% plot corrected velocity
v_north = correct_Deg(:,4);
v_east = correct_Deg(:,5);

plot(num,v_north);
savefig("corrected_vel_north");
plot(num,v_east);
savefig("corrected_vel_east");
% plot(Solutions_NED(:,1), Solutions_NED(:,2))

