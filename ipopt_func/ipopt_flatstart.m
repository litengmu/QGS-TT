clear
clc
tic
warning off

% List of case files (these will be the scripts you want to run)
files = dir('*.m');  % Assuming all cases are *.m files
mpopt(11) = 580;

% Initialize result storage
results_struct = struct('case_name', {}, 'local_optima', {}, 'isglobal', {}, 'isshortest', {}, 'distances', {},...
                        'success_count', {}, 'failure_count', {}, ...
                             'iteration_limit_count', {});  % Add field for average time

% Loop over each case
for file_idx = 1:length(files)
    % Load each case file (this should initialize necessary variables for the case)
    casename = files(file_idx).name(1:end-2);  % Removing ".m" from file name
    fprintf('Running case: %s\n', casename);
    CASE_DEFINE % Assuming each case file defines the required variables (A, B, C, etc.)
    % Define the number of points to be evaluated
    Npnts = 100;

    local_optima = [];
    success_count = 0;
    failure_count = 0;
    iteration_limit_count = 0;
    total_computation_time = 0;
    distances = [];
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
            % Check if the current solution is unique
            new_objective = opf_costfcn(x, om);
            new_x = x;

            % Check uniqueness: Compare both the objective and the variable solution
            is_unique = true;
            for k = 1:length(local_optima)
                if abs(local_optima(k).objective - new_objective) < 1e-4 
                    is_unique = false;  % This solution is a duplicate
                    break;
                end
            end

            % If it's unique, add it to the list of local optima
            if is_unique
                local_optima = [local_optima; struct('objective', new_objective, 'x', new_x)];
            end
        elseif info.status == -1  % Iteration limit reached
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
    results_struct(file_idx).case_name = casename;
    results_struct(file_idx).local_optima = local_optima;
    results_struct(file_idx).isglobal = isglobal;
    results_struct(file_idx).isshortest = isshortest;
    results_struct(file_idx).distances = distances;

    results_struct(file_idx).success_count = success_count;
    results_struct(file_idx).failure_count = failure_count;
    results_struct(file_idx).iteration_limit_count = iteration_limit_count;

    % Display summary for each case
    fprintf('Results for case %s:\n', casename);
    fprintf('Success Count: %d\n', success_count);
    fprintf('Failure Count: %d\n', failure_count);
    fprintf('Iteration Limit Count: %d\n', iteration_limit_count);
    fprintf('Total Iterations: %d\n\n', Npnts);
end

% Save the results to a MAT file for later analysis
save('IPOPT_optimization_results.mat', 'results_struct');

% Display total execution time
toc
