function [Positions, outliers] = cw1_multi_epoch_positioning(Ranges, r_ea_e)
% cw1_multi_epoch_positioning() - compute the position at multiple epochs 
% using the all rows of the pseudo-range measurements. Starting by an 
% initial position at the centre of the Earth.
%
% Reference:
% COMP0130: ROBOT VISION AND NAVIGATION
% Workshop 1: Task 2
%
% This function created 22/02/2021 by Yuhang Zhang
%
% Inputs:
%   Ranges           Pesudo-ranges table read from files.
%   r_ea_e           initial position (in NED).
%
% Outputs:
%   Positions        Computed positions for every epochs with outliers 
%                    removed. (ECEF).
%   outliers         Detected outliers.
%                    Column 1: Epoch
%                    Column 2: iteration
%                    Column 3: satellite number
%                    Column 4: residual
%                    Column 5: index of the measurement that the outlier
%                              is at.

% Copyright 2021, Yuhang Zhang
% License: BSD; see license.txt for details

% =========================================================================
% Begins

% load constants
Define_Constants;

% satellite numbers
numbers = Ranges{1,2:end};
times = Ranges{2:end,1};

prev_offset = 0;
Positions = zeros(4, size(times,1));
outliers = [];


[r_ea_e,v_eb_e] = pv_NED_to_ECEF(r_ea_e(1),r_ea_e(2),r_ea_e(3),[0;0;0]);

for epoch=1:size(times,1)
    disp("Epoch ["+epoch+"]");
    % cartesian ECEF positions of satellites
    r_ej_e = zeros(size(numbers,2), 3);

    for i=1:size(r_ej_e,1)
        [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(times(epoch), numbers(i));
        r_ej_e(i,:) = sat_r_es_e;
    end
    % disp(r_e0_e);

    % initialise line-of-sight unit vector
    u_aj_e = zeros(size(numbers,2), 3);

    % initialise measurement innovation vector
    delta_z = zeros(size(numbers,2), 1);

    % initialise measurement matrix
    H = ones(size(numbers,2), 4);

    x_plus = zeros(4,1);
    clock_offset = prev_offset;

    disp("################ iteration Starts");
    for itr=1:10
        % predict pseudo ranges
        r_aj_minus = zeros(size(numbers,2),1);
        for i=1:size(numbers,2)
            f = eye(3)*r_ej_e(i,:)'-r_ea_e;
            r_a = sqrt(f.'*f);
            for j=1:5
                % Sagnac effect compensation matrix
                C_e_I = [
                    1 omega_ie*r_a/c 0;
                    -omega_ie*r_a/c 1 0;
                    0 0 1
                ];

                f = C_e_I*r_ej_e(i,:)'-r_ea_e;
                r_a = sqrt(f.'*f);
            end
            r_aj_minus(i) = r_a;

            % Sagnac effect compensation matrix
            C_e_I = [
                1 omega_ie*r_aj_minus(i)/c 0;
                -omega_ie*r_aj_minus(i)/c 1 0;
                0 0 1
            ];
            u_a = (C_e_I*r_ej_e(i,:)'-r_ea_e)/r_aj_minus(i);
            u_aj_e(i,:) = u_a';

            delta_z(i) = Ranges{epoch+1,i+1} - r_aj_minus(i) - clock_offset;
            H(i,1:3) = -u_aj_e(i,:);

        end
        
        % =============================================== Outlier Detection
        % residuals vector
        v = (H*inv(H'*H)*H'-eye(8))*delta_z;

        sigma_p = 5;
        % residuals covariance matrix
        C = (eye(8)-H*inv(H'*H)*H') * sigma_p^2;

        th = 6;
        max_residual = 0;
        max_ind = -1;
        idx = [];
        for j=1:8
            if (abs(v(j))>sqrt(C(j,j))*th)
                disp(abs(v(j))+">"+sqrt(C(j,j))*th);
                idx = [idx j];
                
                if (abs(v(j))>max_residual)
                    max_residual = abs(v(j));
                    max_ind = j;
                    
                end
            end
        end
        
        
        if (max_ind > -1)
            outliers = [
                outliers; 
                times(epoch) itr numbers(max_ind) abs(v(max_ind)) max_ind               
            ];
            delta_z(max_ind) = [];
            H(max_ind,:) = [];
        end
        
        

        % predited state vector
        x_minus = [r_ea_e; clock_offset];

        % least squares
        x_plus = x_minus + (H.'*H)\H.' * delta_z;

        err = abs(norm(r_ea_e) - norm(x_plus(1:3)));
        disp("End of iter["+itr+"] error: " + num2str(err,'%.6f'));

        if(err<0.1 || max_ind>-1)
            break;
        end  

        r_ea_e = x_plus(1:3);
    end
    
    Positions(:,epoch) = x_plus;
    
    prev_offset = x_plus(4);
end

% Ends