function [u, v] = evV(beta_s, aCoor_s, evCoor_r, pl_s, q, gamma)
    
    % _s sender
    % _r reciever
    % q must be a column vector
    n = length(evCoor_r); %[2 x n] size
    p = length(beta_s);
    %U_infv = [U_inf; 0]; % Freestream velocity vector
    
    % Initialize matrices and vectors
    usg = zeros(n, p);
    vsg = zeros(n, p);
    uvg = zeros(n, p);
    vvg = zeros(n, p);

    % Velocity at ith evaluation point due to jth panel
    for i = 1:n
        for j = 1:p
            % Rotation matrices
            g2lM = [cos(beta_s(j)) sin(beta_s(j)); -sin(beta_s(j)) cos(beta_s(j))];
            l2gM = g2lM';
    
            % Evaluation point in the local reference frame of the jth panel
            cCoorl = g2lM * (evCoor_r(:, i) - aCoor_s(:, j));
            r1 = norm(cCoorl);
            r2v = cCoorl - [pl_s(j); 0];
            r2 = norm(r2v);
            
            % Calculate angles theta1 and theta2 in the local system
            theta1 = atan2(cCoorl(2), cCoorl(1));
            theta2 = atan2(r2v(2), r2v(1));
            
            % Source/sink influence calculation (local to global)

            Usl = [-(1/(2*pi)) * log(r2/r1); (theta2 - theta1)/(2*pi)];
            
            Usg = l2gM * Usl; % Transform to global coordinates
            usg(i, j) = Usg(1);
            vsg(i, j) = Usg(2);
            
            % Assemble source/sink submatrix
            
            
            % Vortex influence calculation (local to global)
            Uvl = [Usl(2,1); -1 * Usl(1,1)];
            Uvg = l2gM * Uvl; % Transform to global coordinates
            uvg(i, j) = Uvg(1);
            vvg(i, j) = Uvg(2);
        end         
    end
    % Source velocities at n evaluation points
    us = usg * q;
    vs = vsg * q;
    % Vortex velocities at n evaluation points
    uv = sum(uvg, 2) * gamma;
    vv = sum(vvg, 2) * gamma;

    % General velocity at n evaluation points
    u = us + uv;
    v = vs + vv;

end