function airfoilUI()
    % Global airfoil parameters
    fig = uifigure('Name', 'Airfoil Analysis Input', 'Position', [100, 100, 500, 700]);

    % Number of airfoils
    uilabel(fig, 'Text', 'Number of Airfoils:', 'Position', [20, 650, 150, 30]);
    airfoilCountField = uieditfield(fig, 'numeric', 'Position', [200, 650, 100, 30], 'Value', 1);

    % Number of panels
    uilabel(fig, 'Text', 'Number of Panels:', 'Position', [20, 600, 150, 30]);
    panelsField = uieditfield(fig, 'numeric', 'Position', [200, 600, 100, 30], 'Value', 50);

    % Free-stream velocity
    uilabel(fig, 'Text', 'Free-Stream Velocity [m/s]:', 'Position', [20, 550, 150, 30]);
    velocityField = uieditfield(fig, 'numeric', 'Position', [200, 550, 100, 30], 'Value', 10);

    % Ground height
    uilabel(fig, 'Text', 'Ground height [m]:', 'Position', [20, 500, 150, 30]);
    groundField = uieditfield(fig, 'numeric', 'Position', [200, 500, 100, 30], 'Value', 0);

    % Velocity field resolution
    uilabel(fig, 'Text', 'Velocity field resolution:', 'Position', [20, 450, 150, 30]);
    resField = uieditfield(fig, 'numeric', 'Position', [200, 450, 100, 30], 'Value', 200);

    % Generate Fields Button
    uibutton(fig, 'Text', 'Generate Fields', 'Position', [350, 650, 120, 30], ...
        'ButtonPushedFcn', @(btn, event) generateFields(fig, airfoilCountField.Value, panelsField.Value, velocityField.Value, groundField.Value, resField.Value));
end

function generateFields(fig, an, p, U_inf, grh, res)
    % Delete existing dynamic UI components if any
    delete(findall(fig, 'Type', 'uigridlayout'));

    % Create a panel to hold the grid layout
    panel = uipanel(fig, 'Position', [20, 60, 450, 380], 'Title', 'Airfoil Data Input');

    % Adjust grid layout to accommodate airfoil fields in two rows per airfoil
    grid = uigridlayout(panel);
    grid.RowHeight = repmat({30}, 1, an * 2 + 1); % 2 rows per airfoil + 1 header row
    grid.ColumnWidth = {'1x', '1x', '1x', '1x', '1x'}; % 5 equally spaced columns
    grid.Padding = [10 10 10 10];
    grid.RowSpacing = 5;
    grid.ColumnSpacing = 10;

    % Add a label for the section header
    headerLabel = uilabel(grid, 'Text', 'Input Parameters for Each Airfoil', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    headerLabel.Layout.Row = 1;
    headerLabel.Layout.Column = [1 5]; % Spans all columns

    % Fields for each airfoil
    nacaFields = cell(an, 1);
    chordFields = cell(an, 1);
    alphaFields = cell(an, 1);
    xleFields = cell(an, 1);
    yleFields = cell(an, 1);

    % Add fields for each airfoil
    for i = 1:an
        % Row 1: NACA, Chord, AoA
        nacaLabel = uilabel(grid, 'Text', sprintf('Airfoil %d NACA:', i), 'HorizontalAlignment', 'right');
        nacaLabel.Layout.Row = 2 * i;
        nacaLabel.Layout.Column = 1;

        nacaFields{i} = uieditfield(grid, 'text', 'Value', '0012'); % Default: '0012'
        nacaFields{i}.Layout.Row = 2 * i;
        nacaFields{i}.Layout.Column = 2;

        chordLabel = uilabel(grid, 'Text', sprintf('Airfoil %d Chord:', i), 'HorizontalAlignment', 'right');
        chordLabel.Layout.Row = 2 * i;
        chordLabel.Layout.Column = 3;

        chordFields{i} = uieditfield(grid, 'numeric', 'Value', 1.0); % Default: 1.0
        chordFields{i}.Layout.Row = 2 * i;
        chordFields{i}.Layout.Column = 4;

        aoaLabel = uilabel(grid, 'Text', sprintf('Airfoil %d AoA:', i), 'HorizontalAlignment', 'right');
        aoaLabel.Layout.Row = 2 * i;
        aoaLabel.Layout.Column = 5;

        alphaFields{i} = uieditfield(grid, 'numeric', 'Value', 0.0); % Default: 0.0
        alphaFields{i}.Layout.Row = 2 * i;
        alphaFields{i}.Layout.Column = 6;

        % Row 2: Leading Edge X and Y
        xleLabel = uilabel(grid, 'Text', sprintf('Airfoil %d X LE:', i), 'HorizontalAlignment', 'right');
        xleLabel.Layout.Row = 2 * i + 1;
        xleLabel.Layout.Column = 1;

        xleFields{i} = uieditfield(grid, 'numeric', 'Value', 2 * i); % Default: 2 * i
        xleFields{i}.Layout.Row = 2 * i + 1;
        xleFields{i}.Layout.Column = 2;

        yleLabel = uilabel(grid, 'Text', sprintf('Airfoil %d Y LE:', i), 'HorizontalAlignment', 'right');
        yleLabel.Layout.Row = 2 * i + 1;
        yleLabel.Layout.Column = 3;

        yleFields{i} = uieditfield(grid, 'numeric', 'Value', 0.0); % Default: 0.0
        yleFields{i}.Layout.Row = 2 * i + 1;
        yleFields{i}.Layout.Column = 4;
    end
    % Add a submit button outside the grid, below the panel
    submitButton = uibutton(fig, 'Text', 'Submit', 'ButtonPushedFcn', ...
        @(btn, event) processAirfoilData(nacaFields, chordFields, alphaFields, xleFields, yleFields, an, p, U_inf, grh, res));
    submitButton.Position = [200, 20, 100, 30]; % [X, Y, Width, Height]
end

function processAirfoilData(nacaFields, chordFields, alphaFields, xleFields, yleFields, an, p, U_inf, grh, res)
    % Initialize variables to store input data
    naca = cell(an, 1);
    chord = zeros(an, 1);
    alpha = zeros(an, 1);
    xle = zeros(an, 1);
    yle = zeros(an, 1);

    % Initialize variables for airfoil geometry
    aCoor = zeros(2 * an, p + 1);
    cCoor = zeros(2 * an, p);
    nvec = zeros(2 * an, p);
    tvec = zeros(2 * an, p);
    beta = zeros(an, p);
    pl = zeros(an, p);
    gr = 0; % Initialize ground location

    ki = 1:2:2 * an - 1; % Counter for vector matrices

    % Extract user input from UI fields & create airfoil system geometry
    for k = 1:an
        naca{k} = nacaFields{k}.Value; % Airfoil NACA series
        chord(k) = chordFields{k}.Value; % Chord length
        alpha(k) = alphaFields{k}.Value; % Angle of attack
        xle(k) = xleFields{k}.Value; % Leading edge X position
        yle(k) = yleFields{k}.Value; % Leading edge Y position

        % Generate airfoil geometry using nacaS4 function
        [aCoor(ki(k):ki(k)+1, :), cCoor(ki(k):ki(k)+1, :), nvec(ki(k):ki(k)+1, :), ...
            tvec(ki(k):ki(k)+1, :), beta(k, :), pl(k, :)] = nacaS4(naca{k}, chord(k), p, alpha(k), xle(k), yle(k));
    end
    % Mirror case routine
    % Initialize values
    aCoorm = zeros(2*an, p+1);
    betam = zeros(an, p);
    plm = zeros(an, p);
    
    if grh > 0
        % Initialize values
        ylem = zeros(an,1);
        aCoorm = zeros(2*an, p+1);
        betam = zeros(an, p);
        plm = zeros(an, p);
    
        grhr = min(min(aCoor(ki+1,:)));
        gr = grhr - grh;
        for k = 1:an
            ylem(k,1) = gr - abs(yle(k,1) - gr);
            
            [aCoorm(ki(k):ki(k)+1, :), betam(k, :), plm(k, :)] = nacaS4m(naca{k,1}, chord(k,1), p, -1 * alpha(k,1), xle(k,1), ylem(k,1));
        end
    end
    % =============================================================================================================
    % Hess-Smith Solution
    [q, gamma, Gamma, Cl_kj, Cpi, Cl, Cmle, Ui] = HSSolve(alpha, beta, aCoor, cCoor, nvec, tvec, pl, U_inf, an, ki, chord, xle, yle, naca, grh, gr, betam, aCoorm, plm);
    % Create results table
    ResultsTable = table(string(naca), alpha, chord, xle, yle, Cl, Cmle,...
        'VariableNames', {'NACA', 'Alpha', 'Chord', 'LE x coorinate', 'LE y coorinate', 'Cl', 'Cm_LE'});
    disp(ResultsTable);

    % =============================================================================================================
    % Generate mesh for evaluation of velocity
    [evCoor, x_ev, y_ev] = evFluidMesh(cCoor, an, ki, res, grh, gr);
    % Calculate velocity field and plot
    [ug, vg, um] = veloField(an, ki, cCoor, beta, aCoor, evCoor, pl, q, gamma, U_inf, grh, x_ev, y_ev, res, betam, aCoorm, plm);

    % Save variables to the MATLAB workspace
    assignin('base', 'q', q);
    assignin('base', 'gamma', gamma);
    assignin('base', 'Gamma', Gamma);
    assignin('base', 'KuttaJ_Cl', Cl_kj);
    assignin('base', 'Cl', Cl);
    assignin('base', 'Cm_le', Cmle);
    assignin('base', 'Cp', Cpi);
    assignin('base', 'atangentialU', Ui);
    assignin('base', 'u_field', ug);
    assignin('base', 'v_field', vg);
    assignin('base', 'umag_field', um);
    assignin('base', 'U_inf', U_inf);
    assignin('base', 'grh', grh);
    assignin('base', 'naca', naca);
    assignin('base', 'chord', chord);
    assignin('base', 'alpha', alpha);
    assignin('base', 'xle', xle);
    assignin('base', 'yle', yle);
    assignin('base', 'airfoil_coor', aCoor);
    assignin('base', 'controlp_coor', cCoor);
    assignin('base', 'nvec', nvec);
    assignin('base', 'tvec', tvec);
    assignin('base', 'beta', beta);
    assignin('base', 'panel_lenght', pl);

end