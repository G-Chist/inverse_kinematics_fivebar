clear;      % Clears all variables from the workspace
clc;        % Clears the command window
close all;  % Closes all figure windows

% Define parameters
OC = 100;
AP = 120;
PB = 120;
OA = 140;
BC = 110;

Height_Romi = 61;

Box_start_offset = 15;

% Define warehouse parameters
Height_box1 = 182;
Height_box2 = 224;
Height_top = 266;


% Define circle parameters
% r = 50/2;  % Radius of the circle
% x_center = 100;  % X-coordinate of the center
% y_center = 100;  % Y-coordinate of the center

% Create linspace for angle theta from 0 to 2*pi
% theta = linspace(0, 2*pi, 100);  % 100 points around the circle

% Define simulation parameters
Chunks = 2;

% THESE LINSPACES DEFINE THE TRAJECTORY OF THE MECHANISM'S TIP
% An array of linspaces is also a linspace btw
Px_values = [linspace(80, 130, Chunks),...
             linspace(130, 80, Chunks),...
             linspace(80, 80, Chunks),...
             linspace(80, 130, Chunks),...
             linspace(130, 80, Chunks),...
             linspace(80, 80, Chunks),...
             linspace(80, 130, Chunks)];
Py_values = [linspace(Height_box1+Box_start_offset,Height_box1+Box_start_offset, Chunks),...
             linspace(Height_box1+Box_start_offset,Height_box1+Box_start_offset, Chunks),...
             linspace(Height_box1+Box_start_offset,Height_box2+Box_start_offset, Chunks),...
             linspace(Height_box2+Box_start_offset,Height_box2+Box_start_offset, Chunks),...
             linspace(Height_box2+Box_start_offset,Height_box2+Box_start_offset, Chunks),...
             linspace(Height_box2+Box_start_offset,Height_top+Box_start_offset, Chunks),...
             linspace(Height_top+Box_start_offset,Height_top+Box_start_offset, Chunks)];

% Create figure
figure;
hold on;
grid on;
axis equal;
xlim([-100, 200]);  % Set x-axis limits
ylim([0, 300]);    % Set y-axis limits
xlabel('X-axis');
ylabel('Y-axis');
title('Five Bar Mechanism Animation');

% Store path of P
Px_path = [];
Py_path = [];

% Store angles theta2, theta4
Theta2_path = [];
Theta4_path = [];

% Loop over different Px and Py values
for i = 1:length(Px_values)
    Px = Px_values(i);
    Py = Py_values(i);

    % Store the current position of P
    Px_path = [Px_path, Px];
    Py_path = [Py_path, Py];
    
    % Solve for A and B positions
    syms xA yA xB yB;
    
    eqOA = xA^2 + (yA - Height_Romi)^2 == OA^2;
    eqPA = (xA - Px)^2 + (yA - Py)^2 == AP^2;
    eqPB = (xB - Px)^2 + (yB - Py)^2 == PB^2;
    eqBC = (xB - OC)^2 + (yB - Height_Romi)^2 == BC^2;
    
    % Solve the system of equations
    solutionsA = solve([eqOA, eqPA], [xA, yA]);
    solutionsB = solve([eqBC, eqPB], [xB, yB]);
    
    % Convert solutions to numeric values
    solutions_x_a = double(solutionsA.xA);
    solutions_y_a = double(solutionsA.yA);
    solutions_x_b = double(solutionsB.xB);
    solutions_y_b = double(solutionsB.yB);
    
    % Pick correct solutions based on motion direction
    if solutions_x_a(1) > 0
        x_a = solutions_x_a(2);
        y_a = solutions_y_a(2);
    else
        x_a = solutions_x_a(1);
        y_a = solutions_y_a(1);
    end

    if solutions_x_b(1) > OC
        x_b = solutions_x_b(2);
        y_b = solutions_y_b(2);
    else
        x_b = solutions_x_b(1);
        y_b = solutions_y_b(1);
    end
    
    % Define points
    O = [0, Height_Romi];
    A = [x_a, y_a];
    C = [OC, Height_Romi];
    B = [x_b, y_b];
    P = [Px, Py];

    Ox = O(1);
    Oy = O(2);
    Ax = A(1);
    Ay = A(2);
    Bx = B(1);
    By = B(2);
    Cx = C(1);
    Cy = C(2);
    Px = P(1);
    Py = P(2);
 
    % Now calculate angles theta2 and theta4
    theta2 = angle_ABC(C, O, A) % theta2 = angle COA
    theta4 = 180 - angle_ABC(O, C, B) %theta4 = angle BCQ where Q is to the right of point C = 180 degrees - angle OCB

    % Store the current angles
    Theta2_path = [Theta2_path, theta2];
    Theta4_path = [Theta4_path, theta4];

    % Clear previous plot frame
    cla;

    % Plot the path of P
    plot(Px_path, Py_path, 'k-', 'LineWidth', 2, 'Color', 'r');
    hold on;
    
    % Plot points
    plot([O(1), A(1)], [O(2), A(2)], 'k-', 'LineWidth', 2); % Line OA
    plot([B(1), C(1)], [B(2), C(2)], 'k-', 'LineWidth', 2); % Line BC
    plot([O(1), C(1)], [O(2), C(2)], 'k-', 'LineWidth', 2); % Line OC
    plot([P(1), A(1)], [P(2), A(2)], 'k-', 'LineWidth', 2); % Line PA
    plot([P(1), B(1)], [P(2), B(2)], 'k-', 'LineWidth', 2); % Line PB

    % Plot warehouse
    plot([160, 160], [0, Height_box1], 'k-', 'LineWidth', 2, 'Color', 'g');
    plot([130, 160], [Height_box1, Height_box1], 'k-', 'LineWidth', 2, 'Color', 'g');
    plot([160, 160], [0, Height_box2], 'k-', 'LineWidth', 2, 'Color', 'g');
    plot([130, 160], [Height_box2, Height_box2], 'k-', 'LineWidth', 2, 'Color', 'g');
    plot([160, 160], [0, Height_top], 'k-', 'LineWidth', 2, 'Color', 'g');
    plot([130, 160], [Height_top, Height_top], 'k-', 'LineWidth', 2, 'Color', 'g');
    
    % Plot markers
    plot([O(1), A(1), C(1), B(1), P(1)], [O(2), A(2), C(2), B(2), P(2)], 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

    % Label the points
    text(A(1)-25, A(2), ' A', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
    text(B(1)+5, B(2), ' B', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
    text(C(1)+10, C(2), ' C', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
    text(O(1)-25, O(2), ' O', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
    text(P(1)+5, P(2), ' P', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');

    % Define symbolic variables
    syms FAx FAy FBx FBy FOx FOy FCx FCy FPx FPy TO TC

    % Define equations (assuming a force of 5N is applied to the tip in the negative Y direction)
    F = -5;

    eqAPx = FAx + FPx == 0;
    eqAPy = FAy + FPy + F == 0;
    eqAPT = FAx*(Py-Ay) - FAy*(Px-Ax) == 0;

    eqOAx = FAx + FOx == 0;
    eqOAy = FAy + FOy == 0;
    eqOAT = TO - FAy*(Ax-Ox) + FAx*(Ay-Oy) == 0;

    eqBCx = FBx + FCx == 0;
    eqBCy = FBy + FCy == 0;
    eqBCT = TC + FBx*(By-Cy) + FBy*(Cx-Bx) == 0;

    eqSYSx = FOx + FCx == 0;
    eqSYSy = FOy + FCy + F == 0;
    eqSYST = TO + TC + (Px-Ox)*F == 0;

    % Solve system
    vars = [FAx, FAy, FBx, FBy, FCx, FCy, FOx, FOy, FPx, FPy, TO, TC]; 
    [A, b] = equationsToMatrix([eqAPx, eqAPy, eqAPT, eqOAx, eqOAy, eqOAT, eqBCx, eqBCy, eqBCT, eqSYSx, eqSYSy, eqSYST], vars);
    solution = linsolve(A, b);

    % Convert solution to numeric values
    solution_numeric = double(solution);
    
    % Assign values individually
    FAx_value = solution_numeric(1);
    FAy_value = solution_numeric(2);
    FBx_value = solution_numeric(3);
    FBy_value = solution_numeric(4);
    FCx_value = solution_numeric(5);
    FCy_value = solution_numeric(6);
    FOx_value = solution_numeric(7);
    FOy_value = solution_numeric(8);
    FPx_value = solution_numeric(9);
    FPy_value = solution_numeric(10);
    TO_value  = solution_numeric(11);
    TC_value  = solution_numeric(12);
    
    % Display results
    fprintf("FAx = %f\n", FAx_value);
    fprintf("FAy = %f\n", FAy_value);
    fprintf("FBx = %f\n", FBx_value);
    fprintf("FBy = %f\n", FBy_value);
    fprintf("FCx = %f\n", FCx_value);
    fprintf("FCy = %f\n", FCy_value);
    fprintf("FOx = %f\n", FOx_value);
    fprintf("FOy = %f\n", FOy_value);
    fprintf("FPx = %f\n", FPx_value);
    fprintf("FPy = %f\n", FPy_value);
    fprintf("TO  = %f\n", TO_value);
    fprintf("TC  = %f\n", TC_value);


    % Plot the forces
    quiver(Px, Py, 0, F*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Cx, Cy, FCx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Cx, Cy, 0, FCy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Bx, By, FBx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Bx, By, 0, FBy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Ax, Ay, FAx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Ax, Ay, 0, FAy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Ox, Oy, FOx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Ox, Oy, 0, FOy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Px, Py, FPx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Px, Py+F*2, 0, FPy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);

    % Label the torques
    text(Ox-10, Oy-30, ['TO = ' num2str(TO_value)], 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
    % Label the torques
    text(Cx-40, Cy-50, ['TC = ' num2str(TC_value)], 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
        
    % Update frame
    pause(0.0001);
end
hold off;

    % Define points again
    O = [0, Height_Romi];
    A = [x_a, y_a];
    C = [OC, Height_Romi];
    B = [x_b, y_b];
    P = [Px, Py];

    Ox = O(1);
    Oy = O(2);
    Ax = A(1);
    Ay = A(2);
    Bx = B(1);
    By = B(2);
    Cx = C(1);
    Cy = C(2);
    Px = P(1);
    Py = P(2);
    
    % FBD of link OA

    figure;
    hold on;
    grid on;
    axis equal;
    xlim([-100, 200]);  % Set x-axis limits
    ylim([0, 300]);    % Set y-axis limits
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Five Bar Mechanism FBDs');
    % Plot points
    plot([O(1), A(1)], [O(2), A(2)], 'k-', 'LineWidth', 2); % Line OA
    
    % Plot markers
    plot([O(1), A(1)], [O(2), A(2)], 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

    % Label the points
    text(A(1)-25, A(2), ' A', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
    text(O(1)-25, O(2)+10, ' O', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');


    % Plot the forces
    quiver(Ax, Ay, FAx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Ax, Ay, 0, FAy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Ox, Oy, FOx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Ox, Oy, 0, FOy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);

    % Label the torques
    text(Ox-10, Oy-30, ['TO = ' num2str(TO_value)], 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
    hold off;

    % FBD of link AP
    
    figure;
    hold on;
    grid on;
    axis equal;
    xlim([-100, 200]);  % Set x-axis limits
    ylim([0, 300]);    % Set y-axis limits
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Five Bar Mechanism FBDs');
    % Plot points
    plot([A(1), P(1)], [A(2), P(2)], 'k-', 'LineWidth', 2); % Line AP
    
    % Plot markers
    plot([A(1), P(1)], [A(2), P(2)], 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

    % Label the points
    text(A(1)-25, A(2), ' A', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
    text(P(1)-25, P(2)+10, ' P', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');

    % Plot the forces
    quiver(Px, Py, 0, F*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Ax, Ay, FAx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Ax, Ay, 0, FAy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Px, Py, FPx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Px, Py+F*2, 0, FPy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);

    hold off;

    % FBD of link BC
    
    figure;
    hold on;
    grid on;
    axis equal;
    xlim([-100, 200]);  % Set x-axis limits
    ylim([0, 300]);    % Set y-axis limits
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Five Bar Mechanism FBDs');
    % Plot points
    plot([B(1), C(1)], [B(2), C(2)], 'k-', 'LineWidth', 2); % Line BC
    
    % Plot markers
    plot([B(1), C(1)], [B(2), C(2)], 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

    % Label the points
    text(B(1)-25, B(2)+10, ' B', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
    text(C(1)-25, C(2), ' C', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');

    % Plot the forces
    quiver(Cx, Cy, FCx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Cx, Cy, 0, FCy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Bx, By, FBx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Bx, By, 0, FBy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);

    % Label the torques
    text(Cx-40, Cy-50, ['TC = ' num2str(TC_value)], 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');

    hold off;

    % FBD of link OC
    
    figure;
    hold on;
    grid on;
    axis equal;
    xlim([-100, 200]);  % Set x-axis limits
    ylim([0, 300]);    % Set y-axis limits
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Five Bar Mechanism FBDs');
    % Plot points
    plot([O(1), C(1)], [O(2), C(2)], 'k-', 'LineWidth', 2); % Line OC
    
    % Plot markers
    plot([O(1), C(1)], [O(2), C(2)], 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

    % Label the points
    text(O(1)-25, O(2)+10, ' O', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
    text(C(1)-25, C(2), ' C', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');

    % Plot the forces
    quiver(Cx, Cy, FCx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Cx, Cy, 0, FCy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Ox, Oy, FOx_value*2, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(Ox, Oy, 0, FOy_value*2, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);

    % Label the torques
    text(Ox-10, Oy-30, ['TO = ' num2str(TO_value)], 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
    % Label the torques
    text(Cx-40, Cy-50, ['TC = ' num2str(TC_value)], 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');

    hold off;

% Function to calculate angle given three points
function theta_deg = angle_ABC(A, B, C)
    % Extract coordinates
    x1 = A(1); y1 = A(2);
    x2 = B(1); y2 = B(2);
    x3 = C(1); y3 = C(2);

    % Compute vectors AB and BC
    AB = [x1 - x2, y1 - y2];
    BC = [x3 - x2, y3 - y2];

    % Compute dot product and magnitudes
    dot_product = dot(AB, BC);
    mag_AB = norm(AB);
    mag_BC = norm(BC);

    % Compute the angle in radians
    theta_rad = acos(dot_product / (mag_AB * mag_BC));

    % Convert to degrees
    theta_deg = rad2deg(theta_rad);
end
