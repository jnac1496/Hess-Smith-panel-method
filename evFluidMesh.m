
function [evCoor, x_ev, y_ev] = evFluidMesh(cCoor, an, ki, res, grh, gr)
    % Velocity field plot
    % Input ki, an, aCoorp
    
    % Generate mesh for evaluation of velocity
    
    n = length(cCoor);
    kii = 1:n:an*n;
    % Combine all airfoil coordinates into a single array for bounds calculation
    asysCoor = zeros(an*n,2); 
    for k = 1:an
        asysCoor(kii(k):kii(k)+n-1,1) = cCoor(ki(k),:)';
        asysCoor(kii(k):kii(k)+n-1,2) = cCoor(ki(k)+1,:)';
    end
    
    % Define the bounds of the mesh
    x_min = min(asysCoor(:, 1)) - 1; % Add some padding
    x_max = max(asysCoor(:, 1)) + 1;
    if grh > 0
        y_min = gr;
    else
        y_min = min(asysCoor(:, 2)) - 1;
    end
    y_max = max(asysCoor(:, 2)) + 1;
    
    % Generate a mesh grid
    [x, y] = meshgrid(linspace(x_min, x_max, res), linspace(y_min, y_max, res)); % Adjust resolution as needed
    
    % Initialize mask for the valid region
    valid_mask = true(size(x));
    
    % Loop through each airfoil and exclude points inside the airfoil
    for k = 1:an        
        % Use inpolygon to check which mesh points are inside this airfoil
        [in, on] = inpolygon(x, y, cCoor(ki(k), :), cCoor(ki(k)+1, :));
        
        % Exclude the interior points but keep the outline
        valid_mask(in & ~on) = false;
    end
    
    
    % Extract valid points from the mesh
    x_valid = x(valid_mask);
    y_valid = y(valid_mask);
    
    % Combine mesh and outline points
    x_ev = [x_valid; asysCoor(:, 1)];
    y_ev = [y_valid; asysCoor(:, 2)];
    
    evCoor = [x_ev'; y_ev'];
    
    % Visualize the mesh
    figure;
    scatter(x_ev, y_ev, 1, 'filled'); % Adjust marker size for visualization
    axis equal;
    xlabel('x [m]');
    ylabel('y [m]');
    title('Airfoil system Mesh preview');
end