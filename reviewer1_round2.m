%% Replicate Reviewer 1, Comment #3 Analysis
% Date: 2025-07-28

clear; clc; close all;
f = @(x) 2*x.^2 - 16*x + 20; % Objective function
g = @(x) x.^2 - 2*x - 3;     % Inequality constraint function
color_map_choice = winter;
% Recreate Figure 1a (Original Problem)
figure('Name', 'Reviewer 1 - Comment 3 Recreation');
x_range = -2:0.1:5;
hold on;
plot(x_range, f(x_range), 'r', 'LineWidth', 1.5);
plot(x_range, g(x_range), 'b', 'LineWidth', 1.5);
plot(x_range, zeros(size(x_range)), 'k--'); 
feasible_x = -1:0.1:3;
feasible_g_y = g(feasible_x);
%plot(feasible_x, feasible_g_y, 'k', 'LineWidth', 4); 
plot(feasible_x, zeros(size(feasible_x)), 'k', 'LineWidth', 4);
x_opt = 3;
f_opt = f(x_opt);
plot(x_opt, f_opt, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
plot(-1, f(-1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(x_opt, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
plot(-1, 0, 'gx', 'MarkerSize', 10, 'LineWidth', 2);
text(-0.5, 30, 'f(x)', 'Color', 'r');
text(3, 15, 'g(x)', 'Color', 'b');
text(1, 2, 'Feasible Set', 'HorizontalAlignment', 'center');
xlabel('x');
ylabel('Function Value');
title('(a) Original Problem');
grid on;
axis([-2 5 -20 50]);
axis square;
box on;
hold off;



figure (2)
hold on;
[X, SV] = meshgrid(-2:0.1:5, -3:0.1:3);
Z = f(X);
surf(X, SV, Z, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
colormap(winter); 
theta = linspace(0, 2*pi, 200);
x_feasible_3d = 1 + 2 * cos(theta);
sv_feasible_3d = 2 * sin(theta);
z_feasible_3d = f(x_feasible_3d); 
plot3(x_feasible_3d, sv_feasible_3d, z_feasible_3d, 'k', 'LineWidth', 3.5); 
plot3(-1, 0, f(-1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'LineWidth', 2);
plot3(3, 0, f(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'LineWidth', 2);
title('(b) Modified Problem');
xlabel('x');
ylabel('sv');
zlabel('Objective Function f(x)');
grid on;
axis square;
box on;
view(3); 
hold off;

figure (3)
subplot(1, 3, 1);
hold on;
y_vals_for_color = linspace(-20, 50, 100);
[X_color, Y_color] = meshgrid(x_range, y_vals_for_color);
Z_color = f(X_color); 
pcolor(X_color, Y_color, Z_color,'FaceAlpha', 0.6,'EdgeColor', 'none');
shading interp; 
colormap(gca, color_map_choice);
plot(x_range, f(x_range), 'r', 'LineWidth', 1);
plot(x_feasible_3d, zeros(size(z_feasible_3d)), 'k', 'LineWidth', 4);
plot(-1, f(-1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(3, f(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
xlabel('x');
ylabel('f(x)');
title('(c) Projection on x-f(x) Plane');
grid on;
axis([-2 5 -20 50]);
axis square;
box on;
hold off;


subplot(1, 3, 2);
hold on;
surf(X, SV, Z, 'FaceAlpha', 0.6,'EdgeColor', 'none');
view(2); 
plot(x_feasible_3d, sv_feasible_3d, 'k', 'LineWidth', 1);
plot(x_feasible_3d, zeros(size(x_feasible_3d)), 'k', 'LineWidth', 4);
plot(x_opt, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(-1, 0, 'gx', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('x');
ylabel('sv');
title('(d) Projection on x-sv Plane');
grid on;
axis([-2 5 -3 3]);
colormap(gca, color_map_choice);
plot([3 3], [-3 3], 'b--', 'LineWidth', 1);
plot([5 5], [-3 3], 'b--', 'LineWidth', 1);
text(1.5, 2.5, 'Feasible set', 'HorizontalAlignment', 'center');
text(3.1, 2, 'f = -10', 'Color', 'b');
text(5.1, -1.5, 'f = -10', 'Color', 'b');
axis square;
box on;
hold off;


subplot(1, 3, 3);
hold on;
y_vals_for_color = linspace(-20, 50, 100);
[X_color, Y_color] = meshgrid(x_range, y_vals_for_color);
Z_color = f(X_color); 
pcolor(X_color, Y_color, Z_color,'FaceAlpha', 0.6,'EdgeColor', 'none');
shading interp; 
colormap(gca, color_map_choice);
plot(x_range, g(x_range), 'b', 'LineWidth', 1);
plot(x_feasible_3d, zeros(size(z_feasible_3d)), 'k', 'LineWidth', 4);
plot(x_opt, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(-1, 0, 'gx', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('x');
ylabel('g(x)');
title('(f) Projection on x-g(x) Plane');
grid on;
axis([-2 5 -20 50]);
axis square;
box on;
hold off;




%% Reviewer 1, Comment #4a Analysis
% Date: 2025-07-28
clear; clc; close all;
initial_conditions = [
    20,    3;
    20,   -3;
    10,    1;
    10,   -1;
     4,    0.01;
     4,   -0.01;
     1,    0.01;
     1,   -0.01;
     1,    0.0001;
     1,   -0.0001
];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
tspan = [0, 100]; 

% Prepare a table to store results
num_cases = size(initial_conditions, 1);
results = zeros(num_cases, 5); % Columns: x0, sv0, xf, svf, error

fprintf('Replicating Reviewer 1, Comment #4 Simulation Results...\n\n');
fprintf('%8s %8s | %8s %8s | %15s\n', 'x0', 'sv0', 'x_final', 'sv_final', 'Feas. Error');
fprintf('------------------------------------------------------------\n');

for i = 1:num_cases
    y0 = [initial_conditions(i, 1); initial_conditions(i, 2)];
    [~, y_sol] = ode15s(@qgs_system, tspan, y0, options);
    final_point = y_sol(end, :);
    xf = final_point(1);
    svf = final_point(2);
    H_final = xf^2 - 2*xf - 3 + svf^2;
    error_final = abs(H_final);
    results(i, :) = [y0(1), y0(2), xf, svf, error_final];
    fprintf('%8.4f %8.4f | %8.4f %8.4f | %15.4e\n', ...
            results(i,1), results(i,2), results(i,3), results(i,4), results(i,5));
end




%% Reviewer 1, Comment #4b Analysis
% Date: 2025-07-28

num_trials = 5000; 
x_range = [-5, 25]; 
sv_range = [-5, 5]; 
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
tspan = [0, 50]; % Time span for integration

results = zeros(num_trials, 5); % Columns: x0, sv0, xf, svf, error

upper_manifold_count = 0; 
lower_manifold_count = 0; 

fprintf('Starting QGS simulation with %d random initial points...\n', num_trials);


tic; 
for i = 1:num_trials
    x0 = x_range(1) + (x_range(2) - x_range(1)) * rand;
    sv0 = sv_range(1) + (sv_range(2) - sv_range(1)) * rand;
    y0 = [x0; sv0];
    % Solve the ODE system
    [~, y_sol] = ode45(@qgs_system, tspan, y0, options);
    % Get the final point
    final_point = y_sol(end, :);
    xf = final_point(1);
    svf = final_point(2);
    % Calculate the final feasibility error |H(xf, svf)|
    H_final = xf^2 - 2*xf - 3 + svf^2;
    error_final = abs(H_final);
    % Store results
    results(i, :) = [x0, sv0, xf, svf, error_final];
    if svf > 0
        upper_manifold_count = upper_manifold_count + 1;
    else
        lower_manifold_count = lower_manifold_count + 1;
    end
end

total_time = toc; 

fprintf('\n--- Simulation Finished ---\n');
fprintf('Total computation time for %d trials: %.4f seconds\n', num_trials, total_time);
fprintf('Average time per trial: %.6f seconds\n', total_time / num_trials);
fprintf('\nConvergence Statistics:\n');
fprintf('  Converged to upper manifold (sv > 0): %d times (%.1f%%)\n', ...
        upper_manifold_count, 100 * upper_manifold_count / num_trials);
fprintf('  Converged to lower manifold (sv < 0): %d times (%.1f%%)\n', ...
        lower_manifold_count, 100 * lower_manifold_count / num_trials);
fprintf('\nAverage final feasibility error |H(x_f, sv_f)|: %.4e\n', mean(results(:,5)));


% Plot the Final Results
fprintf('Generating plot of final points...\n');
figure(1);
hold on;

theta = linspace(0, 2*pi, 200);
x_circle = 1 + 2 * cos(theta);
sv_circle = 2 * sin(theta);
plot(x_circle, sv_circle, 'r-', 'LineWidth', 2, 'DisplayName', 'Feasible Set (H=0)');

scatter(results(:,3), results(:,4), 30, 'b', 'filled');

axis equal; 
grid on;
title(['Final Positions of ' num2str(num_trials) ' Random Trials']);
xlabel('x');
ylabel('sv');
legend('show', 'Location', 'northwest');
hold off;


%% Definition of the QGS system as a function
% This function defines the differential equations derived by the reviewer
% y(1) = x, y(2) = sv
function dydt = qgs_system(~, y)
    x = y(1);
    sv = y(2);
    H = x^2 - 2*x - 3 + sv^2;
    dxdt = -(2*x - 2) * H;
    dsvdt = -2*sv * H;
    dydt = [dxdt; dsvdt];
end




%%   Reviewer 1 Comment 4
clear; clc; close all;

% Same initial conditions as the reviewer's simulation
initial_conditions = [
    20,    3;
    20,   -3;
    10,    1;
    10,   -1;
     4,    0.01;
     4,   -0.01;
     1,    0.01;
     1,   -0.01;
     1,    0.0001;
     1,   -0.0001
];

% --- IPOPT Problem Setup ---
% Define the function handles for IPOPT
funcs.objective         = @objective;
funcs.gradient          = @gradient;
funcs.constraints       = @constraints;
funcs.jacobian          = @jacobian;
funcs.jacobianstructure = @jacobianstructure;
funcs.hessian           = @hessian;
funcs.hessianstructure  = @hessianstructure;

% Set IPOPT options
options.ipopt.print_level = 0; % Suppress verbose output
options.ipopt.tol         = 1e-9; % Set solver tolerance

% Define the bounds on the constraints. 
% We have one equality constraint: h(x,sv) = 0
options.cl = 0; % Lower bound on constraints
options.cu = 0; % Upper bound on constraints

% --- Run Simulation and Display Results ---
num_cases = size(initial_conditions, 1);
results_ipopt = zeros(num_cases, 5); % Columns: x0, sv0, xf, svf, error

fprintf('Solving Reviewer 1''s Example with IPOPT...\n\n');
fprintf('%8s %8s | %8s %8s | %15s\n', 'x0', 'sv0', 'x_final', 'sv_final', 'Constraint Val.');
fprintf('------------------------------------------------------------\n');

for i = 1:num_cases
    % Set the initial guess for the variables [x; sv]
    x0 = [initial_conditions(i, 1); initial_conditions(i, 2)];
    
    % Call IPOPT
    [x_sol, info] = ipopt(x0, funcs, options);
    
    % Extract results
    xf = x_sol(1);
    svf = x_sol(2);
    
    % The final constraint value, H(xf, svf), should be near zero
    H_final = xf^2 - 2*xf - 3 + svf^2;
    
    % Store and print results
    results_ipopt(i, :) = [x0(1), x0(2), xf, svf, H_final];
    fprintf('%8.4f %8.4f | %8.4f %8.4f | %15.4e\n', ...
            results_ipopt(i,1), results_ipopt(i,2), ...
            results_ipopt(i,3), results_ipopt(i,4), results_ipopt(i,5));
end

% ----------------------------------------------------------------------
% ----- IPOPT Callback Functions -----
% ----------------------------------------------------------------------

function f = objective(z)
    x = z(1);
    f = 0;
end

function g = gradient(z)
    x = z(1);
    g = [0; 0];
end

function c = constraints(z)
    x = z(1);
    sv = z(2);
    c = x^2 - 2*x - 3 + sv^2;
end

function J = jacobianstructure(~)
    J = sparse(ones(1, 2)); 
end

function J = jacobian(z)
    x = z(1);
    sv = z(2);
    % CORRECTION: The Jacobian's VALUE matrix must also be sparse.
    J = sparse([2*x - 2, 2*sv]);
end

function H = hessianstructure(~)
    H = sparse(tril(ones(2,2))); 
end

function H = hessian(~, sigma, lambda)
    Hxx = lambda * (2);
    Hsv_sv = lambda * (2);
    % CORRECTION: The Hessian's VALUE matrix must also be sparse.
    H_dense = [ Hxx,   0;
                  0, Hsv_sv ];
    H = sparse(tril(H_dense));
end