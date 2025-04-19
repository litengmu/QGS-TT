clear
clc
tic
warning off

% Get a list of all case files in the directory
files = dir('*.m');  % Assuming all cases are *.m files

% Pre-allocate variables for storing convergence results
convergence_rates = zeros(length(files), 2); % Column 1 for success, Column 2 for failure

convergence_data = struct('case_name', {}, 'success_count', {}, 'failure_count', {}, ...
                             'iteration_limit_count', {}, 'total_iterations', {}, ...
                             'average_computation_time', {});  % Add field for average time

% Define the number of attempts (iterations) per case
NPS = 1000;

% Loop through each case file
for file_idx = 1:length(files)
    casename = files(file_idx).name(1:end-2);  % Get the name of the case without extension
    fprintf('Running case: %s\n', casename);
    
    % Load the case data (assuming the case files are properly set up)
    mips_opf_define  % Assuming each case file initializes the system data
    
    % Initialize counters for success and failure
    success_count = 0;
    failure_count = 0;
    iteration_limit_count=0;
    total_iterations = 0;
    total_computation_time = 0;

    % Run the optimization process for NPS attempts
    for i = 1:NPS
        % Initialize theta and random VPQ within bounds
        theta = zeros(nb, 1);
        a = uu(nb+1: end) - ll(nb+1: end);
        b = ll(nb+1: end);
        VPQ = rand(nx-nb, 1);
        VPQ = a .* VPQ + b;
        x0 = [theta; VPQ];
       
          % Start timing for the current optimization run
        tic;
        % Run the MIPS solver
        [x, f, success, Output] = mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
        
        % Update total computation time
        total_computation_time = total_computation_time + toc;  % Add elapsed time for this run

        % Update convergence counts based on the success
        total_iterations = total_iterations + 1;
        if success == -1
            failure_count = failure_count + 1;  % Solver failed
        elseif success == 0
            iteration_limit_count = iteration_limit_count + 1;  % Solver reached max iterations
        elseif success == 1
            success_count = success_count + 1;  % Solver converged
        end
    end
    % Calculate average computation time for the current case
    average_computation_time = total_computation_time / NPS;
    % Store convergence results for the current case
    convergence_rates(file_idx, :) = [success_count, failure_count+iteration_limit_count];
    
    % Save detailed results for each case in the structure
    convergence_data(file_idx).case_name = casename;
    convergence_data(file_idx).success_count = success_count;
    convergence_data(file_idx).failure_count = failure_count;
    convergence_data(file_idx).iteration_limit_count = iteration_limit_count;
    convergence_data(file_idx).total_iterations = total_iterations;
    convergence_data(file_idx).average_computation_time = average_computation_time; 

end

% Display the final convergence rates for all cases
disp('Convergence Results for All Cases:');
disp(array2table(convergence_rates, 'VariableNames', {'Success Count', 'Failure Count'}, 'RowNames', {files.name}));

% Save the detailed convergence data (optional)
save('convergence_data.mat', 'convergence_data');

% Display the total execution time
toc
