function [Heading] = cw1_smooth_heading(D)
Define_Constants;
weight = 0.85;
weight = 0.30;
DR = table2array(D);
H = DR(:,7)*deg_to_rad;
R = DR(:,6);
Heading = [H(1)];
heading = H(1);
for i = 2:length(DR)
    heading = weight*H(i)+(1-weight)*(heading+0.5*(R(i-1)+R(i))*0.5);
    Heading = [Heading; heading];
end
Heading = Heading*rad_to_deg;
end