function [aCoor, beta, pl] = nacaS4m(naca, chord, p, alpha, xle, yle)
% nacaS4 returns the coordinates of NACA 4 digits series airfoil, control
% points at mid-lenght of each panel, unit normal and
% tangential vectors of each panel, the angle and lenght of each panel.
%
%   [xu, yu, xl, yl] = nacaS4(naca, chord, n)
%
%   Inputs:
%     naca - naca 4 digits airfoil series STRING
%     chord - airfoil chord lenght
%     p - number of panels to generate naca airfoil
%     alpha - angle of attack
%
%   Output:
%     [arfCoor]' - xy coordinates of the airfoil surface
%     [conCoor]' - xy coordinates of the control points at mid point of each panel
%     [nvec]' - xy coordinates of the unit normal vectors to each panel
%     [tvec]' - xy coordinates of the unit tangential vectors to each panel
%     beta - angle of the panel with respect to the x axis
%     pl - panel lenghts

    % Number of data points for panels
    n = floor(p/2) + 1;
    % Airfoil parameters
    m = str2double(naca(1))/100; % Maximum camber
    p = str2double(naca(2))/10; % Location maximum camber from LE
    t = str2double(naca(3:4))/100; % Maximum thickness
    
    % x points to evaluate airfoil geometry
        % Distribution by semi-circle to increase points locatation near LE and TE
        phi = linspace(pi, 0, n);
        r = 1/2;
        x = r * cos(phi) + r;
    
    % Mean camber line coordinates yc & Thickness distribution yt
    yc = zeros(1, n);
    yt = zeros(1, n);
    theta = zeros(1, n);
    for i = 1:n
        if x(i) <= 1 * p
            if m > 0
            yc(i) = m/p^2 * (2*p*x(i) - x(i)^2);
            theta(i) = atan((2*m)/p^2 * (p - x(i)));
            end
        else
            if m > 0
            yc(i) = m/(1 - p)^2 * ((1 - 2*p) + 2*p*x(i) - x(i)^2);
            theta(i) = atan(2*m/(1 - p)^2 * (p - x(i)));
            end
        end
        
        yt(i) = t/0.2 * (0.2969*sqrt(x(i)) - 0.1260*x(i) - 0.3516*x(i)^2 + 0.2843*x(i)^3 - 0.1036*x(i)^4);
    end
    
    
    % Airfoil coordinates
    xu = chord*(x - yt.*sin(theta));
    yu = - chord*(yc + yt.*cos(theta));
    yu(1, n) = 0;
    
    xl = chord*(x + yt.*sin(theta));
    yl = - chord*(yc - yt.*cos(theta));
    yl(1, n) = 0;
    
    aCoor = [flip(xl) xu(1,2:n) ; flip(yl) yu(1,2:n)];
    
    % Angle of attact airfoil rotation (about LE)
    if alpha ~= 0 
        alpha = alpha*pi/180;
        arfphi = zeros(1, length(aCoor));
        for i = 1:length(aCoor)
            arfphi(i) = atan(aCoor(2,i)/aCoor(1,i));
            aCoor(1,i) = norm(aCoor(:,i))*cos(arfphi(i)-alpha);
            aCoor(2,i) = norm(aCoor(:,i))*sin(arfphi(i)-alpha);
        end
        aCoor(isnan(aCoor)) = [0; 0];
    end

    % Translate airfoil to LE position
    aCoor = aCoor + [xle; yle];
    
    
    % Control points & normal vectors
    
    %cCoor = zeros(2, length(aCoor)-1);
    beta = zeros(1, length(aCoor)-1);
    %nvec = zeros(3, length(aCoor)-1);
    %tvec = zeros(2, length(aCoor)-1);
    pl = zeros(1, length(aCoor)-1);
    for i = 1:length(aCoor)-1
        pv = aCoor(:, i+1) - aCoor(:, i); % Panel vector
        %cCoor(:, i) = pv/2 + aCoor(:, i);
        beta(i) = atan2(pv(2), pv(1));
        pl(i) = norm(pv);
        %tvec(:, i) = pv/pl(i);
        %nvec(:, i) = cross([0;0;1],[tvec(:,i);0]);
    end
    %nvec(3,:) = [];
end
