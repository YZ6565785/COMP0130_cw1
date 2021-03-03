function [x_plus] = cw1_single_epoch_positioning(Ranges)
% cw1_single_epoch_positioning() - compute the position at time 0 using the
% first set of pseudo-range measurements of the pseudo ranges, which is
% processed without using the initial position. Assume an initial position
% at the centre of the Earth.
%
% COMP0130: ROBOT VISION AND NAVIGATION
% Workshop 1: Task 1a
%
% This function created 22/02/2021 by Yuhang Zhang
%
% Inputs:
%   Ranges        Pseudo-ranges read from files
%
% Outputs:
%   x_plus        Cartesian position of body frame w.r.t. ECEF frame

% Copyright 2021, Yuhang Zhang
% License: BSD; see license.txt for details

% =========================================================================
% Begins

% satellite numbers
numbers = Ranges{1,2:end};

% cartesian ECEF positions of satellites
r_e0_e = zeros(size(numbers,2), 3);

for i=1:size(r_e0_e,1)
    [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(0, numbers(i));
    r_e0_e(i,:) = sat_r_es_e;
end
% disp(r_e0_e);

% load constants
Define_Constants;

% initialise predicted user position
r_ea_e = [0;0;0];

% initialise line-of-sight unit vector
u_aj_e = zeros(size(numbers,2), 3);

% initialise measurement innovation vector
delta_z = zeros(size(numbers,2), 1);

% initialise measurement matrix
H = ones(size(numbers,2), 4);

x_plus = zeros(4,1);
clock_offset = 0;

disp("################ iteration Starts");
for itr=1:10
    % predict pseudo ranges
    r_aj_minus = zeros(size(numbers,2),1);
    for i=1:size(numbers,2)
        f = eye(3)*r_e0_e(i,:)'-r_ea_e;
        r_a = sqrt(f.'*f);
        for j=1:5
            % Sagnac effect compensation matrix
            C_e_I = [
                1 omega_ie*r_a/c 0;
                -omega_ie*r_a/c 1 0;
                0 0 1
            ];

            f = C_e_I*r_e0_e(i,:)'-r_ea_e;
            r_a = sqrt(f.'*f);
        end
        r_aj_minus(i) = r_a;

        % Sagnac effect compensation matrix
        C_e_I = [
            1 omega_ie*r_aj_minus(i)/c 0;
            -omega_ie*r_aj_minus(i)/c 1 0;
            0 0 1
        ];
        u_a = (C_e_I*r_e0_e(i,:)'-r_ea_e)/r_aj_minus(i);
        u_aj_e(i,:) = u_a';

        delta_z(i) = Ranges{2,i+1} - r_aj_minus(i) - clock_offset;
        H(i,1:3) = -u_aj_e(i,:);

    end
    
    % predited state vector
    x_minus = [r_ea_e; clock_offset];
    
    % least squares
    x_plus = x_minus + (H.'*H)\H.' * delta_z;
    
    err = abs(norm(r_ea_e) - norm(x_plus(1:3)));
    disp("iter["+itr+"] error: " + num2str(err,'%.6f'));

    if(err<0.1)
        break;
    end  
    
    r_ea_e = x_plus(1:3);
end



% Ends