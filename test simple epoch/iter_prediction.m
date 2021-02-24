% inputs:
%   r_ej_e_list: position vectors for each number of satellites
%   p: pseudo-range table read from file
%   r_true: ground truth position

function [x_hat_plus] = iter_prediction(r_ej_e_list, p, r_true)
    disp("################ Algorithm Starts");
    Define_Constants;

    % initial predicted user position
    r_ea_e = [0;0;0];

    % line-of-sight unit vector
    u_aj_e_list = zeros(8, 3);


    % measurement innovation vector
    delta_z = zeros(8, 1);

    % measurement matrix
    H = ones(8, 4);

    x_hat_plus = zeros(4,1);

    for itr=1:10
        % predited state vector
        x_hat_minus = [r_ea_e; 60];

        % predict pseudo ranges
        r_aj_minus_list = zeros(8,1);
        for itr=1:5
            for i=1:8
            %     disp(":");

                % Sagnac effect compensation matrix
                C_e_I = [
                    1 omega_ie*r_aj_minus_list(i)/c 0;
                    -omega_ie*r_aj_minus_list(i)/c 1 0;
                    0 0 1
                    ];

                f = C_e_I*r_ej_e_list(i,:)'-r_ea_e;
                r_aj_minus_list(i) = sqrt(f.'*f);
            end
        end
        for i=1:8

            % Sagnac effect compensation matrix
            C_e_I = [1 omega_ie*r_aj_minus_list(i)/c 0; -omega_ie*r_aj_minus_list(i)/c 1 0; 0 0 1];
            u_aj_e = (C_e_I*r_ej_e_list(i,:)'-r_ea_e)/r_aj_minus_list(i);
            u_aj_e_list(i,:) = u_aj_e';

            delta_z(i) = p{2,i+1} - r_aj_minus_list(i) - x_hat_minus(4);
            H(i,1:3) = -u_aj_e_list(i,:);

        end
        % least squares
        x_hat_plus = x_hat_minus + (H.'*H)\H.' * delta_z;
        
        r_ea_e = x_hat_plus(1:3);
        

        % disp("r_aj_minus");
        % disp(r_aj_minus_list);
        % 
        % disp("u_aj_e_minus");
        % disp(u_aj_e_list);
        % 
        % disp("x_hat_plus");
        % disp(x_hat_plus);

        err = abs(norm(r_true) - norm(r_ea_e));
        disp("error:");
        disp(err);

        if(err<0.1)
            break;
        end  
    end
end
