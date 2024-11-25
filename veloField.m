function [ug, vg, um] = veloField(an, ki, cCoor, beta, aCoor, evCoor, pl, q, gamma, U_inf, grh, x_ev, y_ev, res, betam, aCoorm, plm)
    

    % Evaluate velocity at each Evaluation point influenced by kth airfoil jth panels
    n = length(pl);
    kii = 1:n:an*n;
    qk = zeros(n,an);
    uk = zeros(length(evCoor),an);
    vk = zeros(length(evCoor),an);
    
    if grh > 0
        ukm = zeros(length(evCoor),an);
        vkm = zeros(length(evCoor),an);
    end
    
    for k = 1:an
        qk(:,k) = q(kii(k):kii(k)+n-1,1);
    end
    for k = 1:an
        [uk(:,k), vk(:,k)] = evV(beta(k,:), aCoor(ki(k):ki(k)+1,:), evCoor, pl(k,:), qk(:,k), gamma(k));
        if grh > 0
            [ukm(:,k), vkm(:,k)] = evV(betam(k,:), aCoorm(ki(k):ki(k)+1,:), evCoor, plm(k,:), qk(:,k), -1 * gamma(k));
        end
    end
    if grh > 0
        ugr = sum(uk,2);
        vgr = sum(vk,2);
        ugm = sum(ukm,2);
        vgm = sum(vkm,2);
    
        ug = ugr + ugm + U_inf;
        vg = vgr + vgm;
        um = sqrt(ug.^2 + vg.^2);
    else
        ug = sum(uk,2) + U_inf;
        vg = sum(vk,2);
        um = sqrt(ug.^2 + vg.^2);
    end
    % Calculate pressure coefficient in the flow field
    cpg = 1 - (um.^2/U_inf^2);
    % ================================================================================================
    % Generate surface plots for velocity field
    % Create a regular grid for interpolation
    x_min = min(x_ev);
    x_max = max(x_ev);
    y_min = min(y_ev);
    y_max = max(y_ev);
    [x_grid, y_grid] = meshgrid(linspace(x_min, x_max, res), linspace(y_min, y_max, res));
    % Interpolate the u,v and velocity magnitud onto the grid
    u_grid = griddata(x_ev, y_ev, ug, x_grid, y_grid, 'natural');
    v_grid = griddata(x_ev, y_ev, vg, x_grid, y_grid, 'natural');
    um_grid = griddata(x_ev, y_ev, um, x_grid, y_grid, 'natural');
    cpg_grid = griddata(x_ev, y_ev, cpg, x_grid, y_grid, 'natural');
    % Mask NaN values (created during interpolation)
    u_grid(isnan(u_grid)) = 0;
    v_grid(isnan(v_grid)) = 0;
    um_grid(isnan(um_grid)) = 0;
    cpg_grid(isnan(cpg_grid)) = 0;
    % Loop through each airfoil and mask points inside the airfoil shapes
    for k = 1:an
        % Find grid points inside the airfoil
        in = inpolygon(x_grid, y_grid, cCoor(ki(k), :), cCoor(ki(k)+1, :));
        % Mask the velocity magnitude inside the airfoil
        u_grid(in) = NaN;
        v_grid(in) = NaN;
        um_grid(in) = NaN;
        cpg_grid(in) = NaN;
    end
    % Visualize the surface plot
    % Create the surface plot
    figure;
    surf(x_grid, y_grid, u_grid, 'EdgeColor', 'none'); % Smooth surface plot
    colormap jet; % Apply colormap
    cb = colorbar; % Add colorbar
    % Add a label to the colorbar
    ylabel(cb, '[m/s]', 'FontSize', 12);
    view(2);
    % Overlay airfoil outlines
    hold on;
    for k = 1:an
        %x_coord = cCoor(ki(k), :);
        %y_coord = cCoor(ki(k)+1, :);
        plot(aCoor(ki(k), :), aCoor(ki(k)+1, :), 'k-', 'LineWidth', 1.5); % Airfoil outlines
    end
    if grh > 0
        plot(linspace(x_min,x_max,50),y_min*ones(1,50), 'k-', 'LineWidth', 2);
    end
    title('Velocity component: u');
    xlabel('X [m]');
    ylabel('Y [m]');
    axis equal;
    hold off

    figure;
    surf(x_grid, y_grid, v_grid, 'EdgeColor', 'none'); % Smooth surface plot
    colormap jet; % Apply colormap
    cb = colorbar; % Add colorbar
    % Add a label to the colorbar
    ylabel(cb, '[m/s]', 'FontSize', 12);
    view(2);
    % Overlay airfoil outlines
    hold on;
    for k = 1:an
        %x_coord = cCoor(ki(k), :);
        %y_coord = cCoor(ki(k)+1, :);
        plot(aCoor(ki(k), :), aCoor(ki(k)+1, :), 'k-', 'LineWidth', 1.5); % Airfoil outlines
    end
    if grh > 0
        plot(linspace(x_min,x_max,50),y_min*ones(1,50), 'k-', 'LineWidth', 2);
    end
    title('Velocity component: v');
    xlabel('X [m]');
    ylabel('Y [m]');
    axis equal;
    hold off

    figure;
    surf(x_grid, y_grid, um_grid, 'EdgeColor', 'none'); % Smooth surface plot
    colormap jet; % Apply colormap
    cb = colorbar; % Add colorbar
    % Add a label to the colorbar
    ylabel(cb, '[m/s]', 'FontSize', 12);
    view(2);
    % Overlay airfoil outlines
    hold on;
    for k = 1:an
        %x_coord = cCoor(ki(k), :);
        %y_coord = cCoor(ki(k)+1, :);
        plot(aCoor(ki(k), :), aCoor(ki(k)+1, :), 'k-', 'LineWidth', 1.5); % Airfoil outlines
    end
    if grh > 0
        plot(linspace(x_min,x_max,50),y_min*ones(1,50), 'k-', 'LineWidth', 2);
    end
    title('Velocity magnitude');
    xlabel('X [m]');
    ylabel('Y [m]');
    axis equal;
    hold off

    figure;
    surf(x_grid, y_grid, cpg_grid, 'EdgeColor', 'none'); % Smooth surface plot
    colormap jet; % Apply colormap
    cb = colorbar; % Add colorbar
    % Add a label to the colorbar
    ylabel(cb, 'Cp', 'FontSize', 12);
    view(2);
    % Overlay airfoil outlines
    hold on;
    for k = 1:an
        %x_coord = cCoor(ki(k), :);
        %y_coord = cCoor(ki(k)+1, :);
        plot(aCoor(ki(k), :), aCoor(ki(k)+1, :), 'k-', 'LineWidth', 1.5); % Airfoil outlines
    end
    if grh > 0
        plot(linspace(x_min,x_max,50),y_min*ones(1,50), 'k-', 'LineWidth', 2);
    end
    title('Pressure coefficient');
    xlabel('X [m]');
    ylabel('Y [m]');
    axis equal;
    hold off    
end