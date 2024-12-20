function [q, gamma, Gamma, Cl_kj, Cpi, Cl, Cmle, Ui] = HSSolve(alpha, beta, aCoor, cCoor, nvec, tvec, pl, U_inf, an, ki, chord, xle, yle, naca, grh, gr, betam, aCoorm, plm) 
    % Assemble linear system
    n = length(cCoor); 
    % Single airfoil cases:
    if an == 1
        % No ground case
        if grh == 0
            [As, av, cs, cv, bs, bv, usg, vsg, uvg, vvg] = HSkk(beta(1,:), aCoor(1:2,:), cCoor(1:2,:), nvec(1:2,:), tvec(1:2,:), pl(1,:), U_inf);
            A = [As av; cs cv];
            b = [bs; bv];
        % Ground case
        else
            [As, av, cs, cv, bs, bv, usg, vsg, uvg, vvg] = HSkk(beta(1,:), aCoor(1:2,:), cCoor(1:2,:), nvec(1:2,:), tvec(1:2,:), pl(1,:), U_inf);
            % [As, av, cs, cv, bs, bv, usg, vsg, uvg, vvg] = HSkl(beta_s, aCoor_s, cCoor_r, nvec_r, tvec_r, pl_s)
            [Asm, avm, csm, cvm, usgm, vsgm, uvgm, vvgm] = HSklm(betam(1,:), aCoorm(1:2,:), cCoor(1:2,:), nvec(1:2,:), tvec(1:2,:), plm(1,:));
            % Assemble A matrix components considering ground effect (mirror airfoil)
            As = As + Asm;
            av = av + avm;
            cs = cs + csm;
            cv = cv + cvm;
            % Source/Sink & Vortex u & v
            usg = usg + usgm;
            vsg = vsg + vsgm;
            uvg = uvg + uvgm;
            vvg = vvg + vvgm;
            % Assemble linear system
            A = [As av; cs cv];
            b = [bs; bv];        
        end
    % Multiple airfoil cases:
    else
        % No ground case
        if grh == 0
            A = zeros(an*n+an);
            b = zeros(an*n+an,1);
            % Source/Sink & Vortex u & v 
            % * Self-influence Diagonal sub-matrices (n x n)
            % * Cross-influence off-Diagonal sub-matrices (n x n)
            % * Rows (size n each) are influenced airfoil & columns (size n each) are the influencer airfoils 
            usg = zeros(an*n);
            vsg = zeros(an*n);
            uvg = zeros(an*n);
            vvg = zeros(an*n);
        
            kii = 1:n:an*n;
            for k = 1:an % For second loop k is the reciver airfoil (rows of A)
                % Self-influece
                [A(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1), A(kii(k):kii(k)+n-1,an*n+k), A(an*n+k,kii(k):kii(k)+n-1), A(an*n+k,an*n+k), b(kii(k):kii(k)+n-1,1), b(an*n+k,1), ...
                    usg(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1), vsg(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1), uvg(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1), vvg(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1)] = ...
                    HSkk(beta(k,:), aCoor(ki(k):ki(k)+1,:), cCoor(ki(k):ki(k)+1, :), nvec(ki(k):ki(k)+1, :), tvec(ki(k):ki(k)+1, :), pl(k,:), U_inf);
                for l = 1:an
                    if l == k % l is the sender airfoil (columns of A)
                        % when sender is sames as reciever, nothing is perfomed
                    else
                        % Cross-influence
                        [A(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1), A(kii(k):kii(k)+n-1,an*n+l), A(an*n+k,kii(l):kii(l)+n-1), A(an*n+k,an*n+l), ...
                            usg(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1), vsg(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1), uvg(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1), vvg(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1)] = ...
                            HSkl(beta(l,:), aCoor(ki(l):ki(l)+1,:), cCoor(ki(k):ki(k)+1,:), nvec(ki(k):ki(k)+1,:), tvec(ki(k):ki(k)+1,:), pl(l,:));
                    end
                end
            end
        % Ground case
        else
            A = zeros(an*n+an);
            Am = zeros(an*n+an);
            b = zeros(an*n+an,1);
            
            % Source/Sink & Vortex u & v:
            % * Self-influence Diagonal sub-matrices of size (n x n)
            % * Cross-influence off-Diagonal sub-matrices (n x n)
            % * Rows (size n each) are influenced airfoil & columns (size n each) are the influencer airfoils 
            % Mirror SS & Vortex: 
            % * Diagonal sub-matrices of size (n x n) are the influence on each Real airfoil by their respective mirror airfoil
            % * Off-diagonal sub-matrices of size (n x n) are the cross influence on each Real airfoil (k) by remaining mirror airfoils (l)
                 
            usg = zeros(an*n);
            usgm = zeros(an*n);
            vsg = zeros(an*n);
            vsgm = zeros(an*n);
            uvg = zeros(an*n);
            uvgm = zeros(an*n);
            vvg = zeros(an*n);
            vvgm = zeros(an*n);
        
            kii = 1:n:an*n;
            for k = 1:an % For second loop k is the reciver airfoil (rows of A)
                % Self-influece (Real airfoil)
                [A(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1), A(kii(k):kii(k)+n-1,an*n+k), A(an*n+k,kii(k):kii(k)+n-1), A(an*n+k,an*n+k), b(kii(k):kii(k)+n-1,1), b(an*n+k,1), ...
                    usg(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1), vsg(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1), uvg(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1), vvg(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1)] = ...
                    HSkk(beta(k,:), aCoor(ki(k):ki(k)+1,:), cCoor(ki(k):ki(k)+1, :), nvec(ki(k):ki(k)+1, :), tvec(ki(k):ki(k)+1, :), pl(k,:), U_inf);
                % Self-influece (Mirror airfoil)
                [Am(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1), Am(kii(k):kii(k)+n-1,an*n+k), Am(an*n+k,kii(k):kii(k)+n-1), Am(an*n+k,an*n+k), ...
                    usgm(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1), vsgm(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1), uvgm(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1), vvgm(kii(k):kii(k)+n-1,kii(k):kii(k)+n-1)] = ...
                    HSkl(betam(k,:), aCoorm(ki(k):ki(k)+1,:), cCoor(ki(k):ki(k)+1, :), nvec(ki(k):ki(k)+1, :), tvec(ki(k):ki(k)+1, :), plm(k,:));
    
                for l = 1:an
                    if l == k % l is the sender airfoil (columns of A)
                        % when sender is sames as reciever, nothing is perfomed
                    else
                        % Cross-influence (Real(index: l) on Real(index: k) airfoils)
                        [A(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1), A(kii(k):kii(k)+n-1,an*n+l), A(an*n+k,kii(l):kii(l)+n-1), A(an*n+k,an*n+l), ...
                            usg(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1), vsg(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1), uvg(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1), vvg(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1)] = ...
                            HSkl(beta(l,:), aCoor(ki(l):ki(l)+1,:), cCoor(ki(k):ki(k)+1,:), nvec(ki(k):ki(k)+1,:), tvec(ki(k):ki(k)+1,:), pl(l,:));
                        % Cross-influence (Mirror(index: l) on Real(index: k) airfoils)
                        [Am(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1), Am(kii(k):kii(k)+n-1,an*n+l), Am(an*n+k,kii(l):kii(l)+n-1), Am(an*n+k,an*n+l), ...
                            usgm(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1), vsgm(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1), uvgm(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1), vvgm(kii(k):kii(k)+n-1,kii(l):kii(l)+n-1)] = ...
                            HSklm(betam(l,:), aCoorm(ki(l):ki(l)+1,:), cCoor(ki(k):ki(k)+1,:), nvec(ki(k):ki(k)+1,:), tvec(ki(k):ki(k)+1,:), plm(l,:));                    
                    end
                end
            end
            % Assemble A and u, v here
            A = A + Am;
            % Source/Sink & Vortex u & v
            usg = usg + usgm;
            vsg = vsg + vsgm;
            uvg = uvg + uvgm;
            vvg = vvg + vvgm;
        end
    end
    % ======================================================================================================================================================================================================
    % Solve linear system
    x = linsolve(A, b);
    % ======================================================================================================================================================================================================
    % Source/Sink (q) and vortex (gamma) intensities, Airfoil total circulation (Gamma) and Kutta Joukowski Theorem Cl

    if an == 1
        q = x(1:end-1);
        gamma = x(end);
        Gamma = sum(pl(1:end))*gamma;
        Cl_kj = 2 * Gamma/(U_inf);
    else
        q = zeros(an*n,1);
        gamma = zeros(an,1);
        Gamma = zeros(an,1);
        Cl_kj = zeros(an,1);
        for k = 1:an
            q(kii(k):kii(k)+n-1,1) = x(kii(k):kii(k)+n-1);
            gamma(k,1) = x(an*n+k);
            Gamma(k,1) = gamma(k)*sum(pl(k:end));
            Cl_kj(k,1) = 2 * Gamma(k,1)/(U_inf);
        end
    end
    % ======================================================================================================================================================================================================
    % Cp, Cl, Cmle Routine
    U_infv = [U_inf; 0];
    if an == 1
        % Calculate velocity induced at each control point due to source strengths
        us_q = usg * q;
        vs_q = vsg * q;
        
        % Calculate induced velocity due to vortex circulation (gamma)
        uv_gamma = gamma * sum(uvg, 2); % Sum across each row (i-th panel)
        vv_gamma = gamma * sum(vvg, 2); % Sum across each row (i-th panel)
        
        
        % Total velocity at each control point (tangential component)
        Ui = zeros(1, n);
        Cpi = zeros(1, n);
        for i = 1:n
            % Calculate the resultant tangential velocity at control point i
            Ui(1, i) = dot(tvec(:, i), U_infv + [us_q(i); vs_q(i)] + [uv_gamma(i); vv_gamma(i)]);
            
            % Pressure coefficient at control point i
            Cpi(1, i) = 1 - (Ui(1, i)^2 / norm(U_infv)^2);
        end
    
        % Airfoils Lift and Momement (LE) Coefficients
        Cl = 0;
        Cmle = 0;
        for i = 1:n
        Cl = Cl + Cpi(1,i) * pl(1,i)/chord * -nvec(2, i);
        Cmle =  Cmle + Cpi(1,i) * pl(1,i)/chord^2 * dot(cross([cCoor(1, i) - xle; cCoor(2, i) - yle; 0], [nvec(:, i); 0]), [0; 0; 1]);
        end
        
        % Identify Stagnation point
        [~, idx] = min(abs(1 - Cpi)); % idx is the index of the stagnation point
        Cpip = [Cpi(1,1:idx), 1, Cpi(1,idx+1:end)];
    
        xstag = (cCoor(1,idx+1) - cCoor(1,idx))/2 + cCoor(1,idx);
        ystag = (cCoor(2,idx+1) - cCoor(2,idx))/2 + cCoor(2,idx);
        cCoorp = [cCoor(:,1:idx), [xstag; ystag], cCoor(:, idx+1:end)];
    
        % Plot Geometry and Cp distribution along x/c
        figure;
        plot(aCoor(1, :), aCoor(2, :), 'DisplayName', append('NACA: ',naca{1}, ' \alpha = ', sprintf('%.2f', alpha)));
        hold on
        plot(xstag, ystag, 'o', 'DisplayName', 'Stagnation point');
        legend('show', 'Location', 'southoutside');
        xlabel('x [m]');
        ylabel('y [m]');
        if grh > 0
            plot(aCoorm(1, :), aCoorm(2, :));
            plot(linspace(xle,max(aCoor(1,:)),50),gr*ones(1,50),'DisplayName', 'Ground', "Color", 'k', 'LineWidth', 2);
        end
        axis equal
        hold off

        figure;
        plot(cCoorp(1,n/2+2:end)./chord, Cpip(n/2+2:end), 'DisplayName', 'Upper Surface');
        hold on
        plot(cCoorp(1,1:n/2+1)./chord, Cpip(1:n/2+1), 'DisplayName', 'Lower Surface');
        set(gca, 'YDir', 'reverse'); % Invert y-axis for conventional Cp plotting
        xlabel('x / chord');
        ylabel('C_p');
        title('Pressure Coefficient Distribution (C_p)');
        legend('show', 'Location', 'southoutside');
        grid on;
        hold off;
    else
        kii = 1:n:an*n;
        % Sum all velocity contributions at each control point of each airfoil
        us_t = zeros(height(usg), length(usg)/an);
        vs_t = zeros(height(vsg), length(vsg)/an);
        uv_t = zeros(height(uvg), length(uvg)/an);
        vv_t = zeros(height(vvg), length(vvg)/an);    
        for k = 1:an
            us_t = us_t + usg(:,kii(k):kii(k)+n-1);
            vs_t = vs_t + vsg(:,kii(k):kii(k)+n-1);
            uv_t = uv_t + uvg(:,kii(k):kii(k)+n-1);
            vv_t = vv_t + vvg(:,kii(k):kii(k)+n-1);
        end
    
        % Calculate velocity induced at each control point due to source strengths
        us_q = zeros(height(us_t),1);
        vs_q = zeros(height(vs_t),1);
        for k = 1:an
            us_q(kii(k):kii(k)+n-1,1) = us_t(kii(k):kii(k)+n-1,:) * q(kii(k):kii(k)+n-1);
            vs_q(kii(k):kii(k)+n-1,1) = vs_t(kii(k):kii(k)+n-1,:) * q(kii(k):kii(k)+n-1);
        end
    
        % Calculate induced velocity due to vortex circulation (gamma)
        uv_gamma = sum(uv_t, 2);
        vv_gamma = sum(vv_t, 2);    
        for k = 1:an
            uv_gamma(kii(k):kii(k)+n-1,1) = gamma(k) * uv_gamma(kii(k):kii(k)+n-1,1);
            vv_gamma(kii(k):kii(k)+n-1,1) = gamma(k) * vv_gamma(kii(k):kii(k)+n-1,1);
        end
        
        % Total velocity at each control point (tangential component)
        Ui = zeros(an, n);
        Cpi = zeros(an, n);
        for k = 1:an
            for i = 1:n
                % Calculate the resultant tangential velocity at control point i of airfoil k
                Ui(k, i) = dot(tvec(ki(k):ki(k)+1, i), U_infv + [us_q(i); vs_q(i)] + [uv_gamma(i); vv_gamma(i)]);
                
                % Pressure coefficient at control point i
                Cpi(k, i) = 1 - (Ui(k, i)^2 / norm(U_infv)^2);
            end
        end
    
        % Airfoils Lift and Momement (LE) Coefficients
        Cl = zeros(an,1);
        Cmle = zeros(an,1);
        idx = zeros(an,1);
        % For plotting and estimating stagnation point location
        Cpip = zeros(an, length(Cpi)+1); 
        stgCoor = zeros(2*an,1);
        cCoorp = zeros(height(cCoor), length(cCoor)+1);
        for k = 1:an
            for i = 1:n
            Cl(k,1) = Cl(k,1) + Cpi(k,i) * pl(k,i)/chord(k) * -nvec(ki(k)+1, i);
            Cmle(k,1) =  Cmle(k,1) + Cpi(k,i) * pl(k,i)/chord(k)^2 * dot(cross([cCoor(ki(k), i) - xle(k,1); cCoor(ki(k)+1, i) - yle(k,1); 0], [nvec(ki(k):ki(k)+1, i); 0]), [0; 0; 1]);
            end
            % Identify Stagnation point
            [~, idx(k,1)] = min(abs(1 - Cpi(k,:))); % idx is the index of the stagnation point
            Cpip(k,:) = [Cpi(k, 1:idx(k)), 1, Cpi(k,idx(k)+1:end)];
            stgCoor(ki(k):ki(k)+1,1) = (cCoor(ki(k):ki(k)+1,idx(k)+1) - cCoor(ki(k):ki(k)+1,idx(k)))./2 + cCoor(ki(k):ki(k)+1,idx(k));
            cCoorp(ki(k):ki(k)+1, :) = [cCoor(ki(k):ki(k)+1, 1:idx(k)) stgCoor(ki(k):ki(k)+1,1) cCoor(ki(k):ki(k)+1, idx(k)+1:end)];
        end
        % Plot Geometry
        figure;
        hold on
        for k = 1:an
            plot(aCoor(ki(k), :), aCoor(ki(k)+1, :), 'DisplayName', append(sprintf('Airfoil %d: ', k), 'NACA ',naca{k}, ' \alpha = ', sprintf('%.2f', alpha(k))));
            plot(stgCoor(ki(k),1),stgCoor(ki(k)+1,1), 'o', 'DisplayName', append(sprintf('Airfoil %d ', k), 'stagnation point'));
        end
        % Mirror and ground plot
        if grh > 0
            for k = 1:an
                plot(aCoorm(ki(k), :), aCoorm(ki(k)+1, :),'DisplayName', sprintf('Mirror Airfoil %d ', k),'Color', [.5 .5 .5]);
            end
        plot(linspace(min(min(aCoor(ki,:))),max(max(aCoor(ki,:))),50),gr*ones(1,50),'DisplayName', 'Ground', "Color", 'k', 'LineWidth', 2);
        end
        legend('show', 'Location', 'southoutside');
        xlabel('x [m]');
        ylabel('y [m]');
        title('Airoils system');
        axis equal;
        hold off
        % Plot Cp distribution along x/c
        colors = lines(an);
        figure;
        hold on
        for k = 1:an
            plot((cCoorp(ki(k),n/2+2:end) - xle(k))./chord(k), Cpip(k,n/2+2:end), 'DisplayName', sprintf('Upper Surface - Airfoil %d', k), 'LineStyle','--', 'Color', colors(k,:));
            plot((cCoorp(ki(k),1:n/2+1) - xle(k))./chord(k), Cpip(k,1:n/2+1), 'DisplayName', sprintf('Lower Surface - Airfoil %d', k), 'LineStyle','-', 'Color', colors(k,:));
            set(gca, 'YDir', 'reverse'); % Invert y-axis for conventional Cp plotting
        end
        xlabel('x / chord');
        ylabel('C_p');
        title('Pressure Coefficient Distribution (C_p)');
        legend('show', 'Location', 'southoutside');
        grid on;
    end
end
