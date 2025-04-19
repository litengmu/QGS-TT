clear
clc
tic
warning off

% List of case files (these will be the scripts you want to run)
files = dir('*.m');  % Assuming all cases are *.m files
mpopt(11)=580;

% Initialize result storage
convergence_results = struct('case_name', {}, 'success_count', {}, 'failure_count', {}, ...
                             'iteration_limit_count', {}, 'total_iterations', {}, ...
                             'average_computation_time', {});  % Add field for average time


% Loop over each case
for file_idx = 1:length(files)
    % Load each case file (this should initialize necessary variables for the case)
    casename = files(file_idx).name(1:end-2);  % Removing ".m" from file name
    fprintf('Running case: %s\n', casename);
    CASE_DEFINE % Assuming each case file defines the required variables (A, B, C, etc.)

    % Initialize counters for success and failure
    success_count = 0;
    failure_count = 0;
    iteration_limit_count = 0;
    total_iterations = 0;
     total_computation_time = 0;
    % Define the number of points to be evaluated
    Npnts = 1000;

    xv=[];
    success_rate=0;
    nA = size(A, 1);                %% number of original linear constraints
    nx = length(x0);
    f = branch(:, F_BUS);                           %% list of "from" buses
    t = branch(:, T_BUS);                           %% list of "to" buses
    Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
    Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
    Cl = Cf + Ct;
    Cb = Cl' * Cl + speye(nb);
    Cl2 = Cl(il, :);
    Cg = sparse(gen(:, GEN_BUS), (1:ng)', 1, nb, ng);
    nz = nx - 2*(nb+ng);
    nxtra = nx - 2*nb;
    Js = [
        Cb      Cb      Cg              sparse(nb,ng)   sparse(nb,nz);
        Cb      Cb      sparse(nb,ng)   Cg              sparse(nb,nz);
        Cl2     Cl2     sparse(nl2, 2*ng)               sparse(nl2,nz);
        Cl2     Cl2     sparse(nl2, 2*ng)               sparse(nl2,nz);
        A;
        ];
    [f, df, d2f] = opf_costfcn(x0, om);
    Hs = tril(d2f + [
        Cb  Cb  sparse(nb,nxtra);
        Cb  Cb  sparse(nb,nxtra);
        sparse(nxtra,nx);
        ]);

    %% set options struct for IPOPT
    options.ipopt = ipopt_options([], mpopt);

    %% extra data to pass to functions
    options.auxdata = struct( ...
        'om',       om, ...
        'Ybus',     Ybus, ...
        'Yf',       Yf(il,:), ...
        'Yt',       Yt(il,:), ...
        'mpopt',    mpopt, ...
        'il',       il, ...
        'A',        A, ...
        'nA',       nA, ...
        'neqnln',   2*nb, ...
        'niqnln',   2*nl2, ...
        'Js',       Js, ...
        'Hs',       Hs    );

    %% define variable and constraint bounds
    options.lb = ll;
    options.ub = uu;
    options.cl = [zeros(2*nb, 1); -Inf*ones(2*nl2, 1); l];
    options.cu = [zeros(2*nb, 1);     zeros(2*nl2, 1); u];

    %% assign function handles
    funcs.objective         = @objective;
    funcs.gradient          = @gradient;
    funcs.constraints       = @constraints;
    funcs.jacobian          = @jacobian;
    funcs.hessian           = @hessian;
    funcs.jacobianstructure = @(d) Js;
    funcs.hessianstructure  = @(d) Hs;
    %funcs.jacobianstructure = @jacobianstructure;
    %funcs.hessianstructure  = @hessianstructure;

    % Loop over the points for optimization
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

        % Run IPOPT optimization
        [x, info] = ipopt(y0, funcs, options);

        % Update total computation time
        total_computation_time = total_computation_time + toc;  % Add elapsed time for this run

        % Check convergence status and update counters
        if info.status == 0 || info.status == 1  % Converged successfully
            success_count = success_count + 1;
          
        elseif info.status == -1  % Iteration limit reached
            iteration_limit_count = iteration_limit_count + 1;
        else  % Computational failure
            failure_count = failure_count + 1;
        end

        total_iterations = total_iterations + 1;  % Sum of all iterations for this run

    end
    % Calculate average computation time for the current case
    average_computation_time = total_computation_time / Npnts;
    % Store the results for the current case
    convergence_results(file_idx).case_name = casename;
    convergence_results(file_idx).success_count = success_count;
    convergence_results(file_idx).failure_count = failure_count;
    convergence_results(file_idx).iteration_limit_count = iteration_limit_count;
    convergence_results(file_idx).total_iterations = total_iterations;
    convergence_results(file_idx).average_computation_time = average_computation_time; 

    % Display summary for each case
    fprintf('Results for case %s:\n', casename);
    fprintf('Success Count: %d\n', success_count);
    fprintf('Failure Count: %d\n', failure_count);
    fprintf('Iteration Limit Count: %d\n', iteration_limit_count);
    fprintf('Total Iterations: %d\n\n', total_iterations);
end

% Save the results to a MAT file for later analysis
save('convergence_results.mat', 'convergence_results');

% Display total execution time
toc


