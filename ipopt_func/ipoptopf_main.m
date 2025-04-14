
clc
clear global
clear

warning off
tic
CASE_DEFINE
mpopt(11)=580;


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

% %% check Jacobian and Hessian structure
% xr                  = rand(size(x0));
% lambda              = rand(2*nb+2*nl2, 1);
% options.auxdata.Js  = jacobian(xr, options.auxdata);
% options.auxdata.Hs  = tril(hessian(xr, 1, lambda, options.auxdata));
% Js1 = options.auxdata.Js;
% options.auxdata.Js = Js;
% Hs1 = options.auxdata.Hs;
% [i1, j1, s] = find(Js);
% [i2, j2, s] = find(Js1);
% if length(i1) ~= length(i2) || norm(i1-i2) ~= 0 || norm(j1-j2) ~= 0
%     error('something''s wrong with the Jacobian structure');
% end
% [i1, j1, s] = find(Hs);
% [i2, j2, s] = find(Hs1);
% if length(i1) ~= length(i2) || norm(i1-i2) ~= 0 || norm(j1-j2) ~= 0
%     error('something''s wrong with the Hessian structure');
% end

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


Npnts=10000;
%%-----  run opf  -----
for j=1:Npnts
    theta = zeros(nb, 1);
    a = uu(nb+1: nx) - ll(nb+1: nx);   %%上下限间距
    b = ll(nb+1: nx);                   %%下限
    VPQ = rand(nx-nb, 1);                % VPQ = lhsdesign(nx-nb, 1);
    VPQ = a.* VPQ + b;                   %%（uu-ll）*(0~1)+ll?   分层随机抽点！
    y0(:,j) = [theta; VPQ];

    %% run the optimization
    [x, info] = ipopt(y0(:,j),funcs,options);

    if info.status == 0 || info.status == 1
        success = 1;
        success_rate=success_rate+1;
        xv=[xv x];
    else
        success = 0;
    end
    if isfield(info, 'iter')
        output.iterations = info.iter;
    else
        output.iterations = [];
    end
end

toc
success_rate=success_rate/Npnts;




