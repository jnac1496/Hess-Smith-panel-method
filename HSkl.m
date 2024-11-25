function [As, av, cs, cv, usg, vsg, uvg, vvg] = HSkl(beta_s, aCoor_s, cCoor_r, nvec_r, tvec_r, pl_s)
    
    % _s sender
    % _r reciever 
    n = length(cCoor_r);
    %U_infv = [U_inf; 0]; % Freestream velocity vector
    
    % Initialize matrices and vectors
    usg = zeros(n);
    vsg = zeros(n);
    uvg = zeros(n);
    vvg = zeros(n);
    As = zeros(n);
    Av = zeros(n);
    %bs = zeros(n, 1);
    cs = zeros(1, n);
    cvv = zeros(1, n);
    av = zeros(n, 1);
    % No flow through walls boundary condition
    for i = 1:n
        for j = 1:n
            % Rotation matrices
            g2lM = [cos(beta_s(j)) sin(beta_s(j)); -sin(beta_s(j)) cos(beta_s(j))];
            l2gM = g2lM';
    
            % Control point in the local reference frame of the jth panel
            cCoorl = g2lM * (cCoor_r(:, i) - aCoor_s(:, j));
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
            As(i, j) = dot(nvec_r(:, i), [usg(i, j); vsg(i, j)]);
            
            % Vortex influence calculation (local to global)
            Uvl = [Usl(2,1); -1 * Usl(1,1)];
            Uvg = l2gM * Uvl; % Transform to global coordinates
            uvg(i, j) = Uvg(1);
            vvg(i, j) = Uvg(2);
            
            % Assemble vortex submatrix
            Av(i, j) = dot(nvec_r(:, i), [uvg(i, j); vvg(i, j)]);
        end
        
        % Assemble vortex vector av
        av(i) = sum(Av(i, :));
    end
    
    % Implement the Kutta condition (sum contributions at the trailing edge)
    for j = 1:n
        cs(j) = dot(tvec_r(:, 1), [usg(1, j); vsg(1, j)]) + dot(tvec_r(:, end), [usg(end, j); vsg(end, j)]);
        cvv(j) = dot(tvec_r(:, 1), [uvg(1, j); vvg(1, j)]) + dot(tvec_r(:, end), [uvg(end, j); vvg(end, j)]);
    end
    
    % Calculate vortex scalar cv
    cv = sum(cvv);
    
end