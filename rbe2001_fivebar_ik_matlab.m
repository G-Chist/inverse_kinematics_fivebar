% Define parameters
OC = 40;
AP = 50;
PB = 50;
OA = 60;
BC = 60;

% Define Px and Py range for animation
Px_values = linspace(-20, 70, 100); % Varying Px
Py_values = linspace(60, 60, 100); % Varying Py

% Create figure
figure;
hold on;
grid on;
axis equal;
xlim([-100, 200]);  % Set x-axis limits
ylim([0, 200]);    % Set y-axis limits
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
    
    eqOA = xA^2 + yA^2 == OA^2;
    eqPA = (xA - Px)^2 + (yA - Py)^2 == AP^2;
    eqPB = (xB - Px)^2 + (yB - Py)^2 == PB^2;
    eqBC = (xB - OC)^2 + yB^2 == BC^2;
    
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

    if solutions_x_b(1) < OC
        x_b = solutions_x_b(2);
        y_b = solutions_y_b(2);
    else
        x_b = solutions_x_b(1);
        y_b = solutions_y_b(1);
    end
    
    % Define points
    O = [0, 0];
    A = [x_a, y_a];
    C = [OC, 0];
    B = [x_b, y_b];
    P = [Px, Py];

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
    
    % Plot markers
    plot([O(1), A(1), C(1), B(1), P(1)], [O(2), A(2), C(2), B(2), P(2)], 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

    % Update frame
    pause(0.01);
end

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