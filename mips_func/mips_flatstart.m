clear
clc
tic
warning off

% List of case files (these will be the scripts you want to run)
files = dir('*.m');  % Assuming all cases are *.m files
mpopt(11) = 560;

% Initialize result storage
MIPS_results_struct = struct('case_name', {}, 'local_optima', {}, 'isglobal', {}, 'isshortest', {}, 'distances', {},...
                        'success_count', {}, 'failure_count', {}, ...
                         'iteration_limit_count', {},'distances_controlled',{});  % Add field for average time

% Loop over each case
for file_idx = 1:length(files)
    % Load each case file (this should initialize necessary variables for the case)
    casename = files(file_idx).name(1:end-2);  % Removing ".m" from file name
    fprintf('Running case: %s\n', casename);
    mips_opf_define % Assuming eac
    % h case file defines the required variables (A, B, C, etc.)
    % Define the number of points to be evaluated
    Npnts = 1000;

    local_optima = [];
    success_count = 0;
    failure_count = 0;
    iteration_limit_count = 0;
    total_computation_time = 0;
    distances = [];

    for j = 1:Npnts
        % Initialize random point (for each run)
        theta = zeros(nb, 1);
        a = uu(nb+1: nx) - ll(nb+1: nx);   % upper - lower bounds
        b = ll(nb+1: nx);                   % lower bounds
        VPQ = rand(nx-nb, 1);
        VPQ = a .* VPQ + b;  % Random sampling within bounds
        y0 = [theta; VPQ];

        % Start timing for the current optimization run
        tic;

    [x, f, success, Output] = ...
        mips(f_fcn, y0, A, l, u, xmin, xmax, gh_fcn, hess_fcn,opt);
  
        % Update total computation time
        total_computation_time = total_computation_time + toc;  % Add elapsed time for this run

        % Check convergence status and update counters
        if success == 1  % Converged successfully
            success_count = success_count + 1;
            % Check if the current solution is unique
            new_x = x;
            % Check uniqueness: Compare both the objective and the variable solution
            is_unique = true;
            for k = 1:length(local_optima)
                if abs(local_optima(k).objective - f)/f < 1e-3 
                    is_unique = false;  % This solution is a duplicate
                    break;
                end
            end

            % If it's unique, add it to the list of local optima
            if is_unique
                local_optima = [local_optima; struct('objective', f, 'x', new_x)];
            end
        elseif success == 0  % Iteration limit reached
            iteration_limit_count = iteration_limit_count + 1;
        else  % Computational failure
            failure_count = failure_count + 1;
        end
    end

    % After running all 1000 points, get the best objective (global optimum)
    best_objective = min([local_optima.objective]);
    best_x = local_optima([local_optima.objective] == best_objective).x;

    % Run the flat-start point OPF for comparison
    flat_start_results = runopf(casename);
    flat_start_objective = flat_start_results.f;
    flat_start_x = flat_start_results.x0;

    % Check if the flat-start point is the global optimum
    isglobal = flat_start_objective -best_objective<1e-4;

    % Compute distance between the flat-start point and each local optimum
    for k = 1:size(local_optima,1)
        distance = norm(flat_start_x - local_optima(k).x);
        distances = [distances; distance];
    end

    % Check if the flat-start point is the closest to the global optimum
    min_distance = min(distances);
    isshortest = (min_distance == norm(flat_start_x - best_x));

    % Store the results for the current case
    MIPS_results_struct(file_idx).case_name = casename;
    MIPS_results_struct(file_idx).local_optima = local_optima;
    MIPS_results_struct(file_idx).isglobal = isglobal;
    MIPS_results_struct(file_idx).isshortest = isshortest;
    MIPS_results_struct(file_idx).distances = distances;

    MIPS_results_struct(file_idx).success_count = success_count;
    MIPS_results_struct(file_idx).failure_count = failure_count;
    MIPS_results_struct(file_idx).iteration_limit_count = iteration_limit_count;

    % Display summary for each case
    fprintf('Results for case %s:\n', casename);
    fprintf('Success Count: %d\n', success_count);
    fprintf('Failure Count: %d\n', failure_count);
    fprintf('Iteration Limit Count: %d\n', iteration_limit_count);
    fprintf('Total Iterations: %d\n\n', Npnts);
end

% Save the results to a MAT file for later analysis
save('MIPS_optimization_results.mat', 'MIPS_results_struct');

% Display total execution time
toc

%% 仅仅看控制变量的结果

for file_idx = 1:length(files)
    distances_controlled=[];
    % Load each case file (this should initialize necessary variables for the case)
    casename = files(file_idx).name(1:end-2);  % Removing ".m" from file name
    fprintf('Running case: %s\n', casename);
    mips_opf_define % Assuming eac

    for k = 1:size(MIPS_results_struct(file_idx).local_optima,1)
        distance_controlled = norm(MIPS_results_struct(file_idx).local_optima(k).x(nb+1:end) - x0(nb+1:end));
        distances_controlled = [distances_controlled; distance_controlled];
    end
    MIPS_results_struct(file_idx).distances_controlled = distances_controlled;
end






