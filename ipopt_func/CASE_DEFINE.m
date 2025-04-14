
casename='WB2';
casedata =feval(casename);

casedata.branch(casedata.branch(:,6)<1e-4,6) = 9999;  
% fprintf('-------Start Solve Case: %s....\n\n',casename)

%% define named indices into data matrices

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;



%% add zero columns to bus, gen, branch for multipliers, etc if needed
nb   = size(casedata.bus, 1);    %% number of buses
nl   = size(casedata.branch, 1); %% number of branches
ng   = size(casedata.gen, 1);    %% number of dispatchable injections
if size(casedata.bus,2) < MU_VMIN
    casedata.bus = [casedata.bus zeros(nb, MU_VMIN-size(casedata.bus,2)) ];
end
if size(casedata.gen,2) < MU_QMIN
    casedata.gen = [ casedata.gen zeros(ng, MU_QMIN-size(casedata.gen,2)) ];
end
if size(casedata.branch,2) < MU_ANGMAX
    casedata.branch = [ casedata.branch zeros(nl, MU_ANGMAX-size(casedata.branch,2)) ];
end

mpopt = mpoption('OPF_ALG',580,'OUT_ALL',0); 

mpc = ext2int(casedata);


om = opf_setup(mpc, mpopt); 
om = build_cost_params(om);  

%% ----- initialization -----

verbose = mpopt(31);    
feastol = mpopt(81);    %% PDIPM_FEASTO
gradtol = mpopt(82);    %% PDIPM_GRADTOL
comptol = mpopt(83);    %% PDIPM_COMPTOL
costtol = mpopt(84);    %% PDIPM_COSTTOL
max_it  = mpopt(85);    %% PDIPM_MAX_IT
max_red = mpopt(86);    %% SCPDIPM_RED_IT
step_control = (mpopt(11) == 565);  %% OPF_ALG == 565, MIPS-sc
%%565 { MIPS-sc, step-controlled variant of MIPS 

if feastol == 0
    feastol = mpopt(16);    %% = OPF_VIOLATION by default，5e-6
end

opt = struct(   'feastol', feastol, ...
                'gradtol', gradtol, ...
                'comptol', comptol, ...
                'costtol', costtol, ...
                'max_it', 200, ...
                'max_red', 1, ...
                'step_control', step_control, ...
                'cost_mult', 1e-4, ...
                'verbose', verbose  );

%% unpack data
mpc = get_mpc(om);  %om.mpc
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
[vv, ll, nn] = get_idx(om);
%% problem dimensions 
global nb nl ng nl2
nb = size(bus, 1);          
nl = size(branch, 1);       
ng = size(gen, 1);    
ny = getN(om, 'var', 'y');  
nx = 2*nb + 2*ng;    

%% linear constraints
 [A, l, u] = linear_constraints(om);   

%% bounds on optimization vars
[x0, xmin, xmax] = getv(om);  

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% try to select an interior initial point  
ll = xmin; uu = xmax;
ll(xmin == -Inf) = -1e10;   
uu(xmax ==  Inf) =  1e10;
x0 = (ll + uu) / 2;  
Varefs = bus(bus(:, BUS_TYPE) == REF, VA) * (pi/180);  
x0(vv.i1.Va:vv.iN.Va) = Varefs(1);  

if ny > 0
    ipwl = find(gencost(:, MODEL) == PW_LINEAR);
%     PQ = [gen(:, PMAX); gen(:, QMAX)];
%     c = totcost(gencost(ipwl, :), PQ(ipwl));
    c = gencost(sub2ind(size(gencost), ipwl, NCOST+2*gencost(ipwl, NCOST)));    %% largest y-value in CCV data
    x0(vv.i1.y:vv.iN.y) = max(c) + 0.1 * abs(max(c));
%     x0(vv.i1.y:vv.iN.y) = c + 0.1 * abs(c);
end

%% find branches with flow limits 
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
nl2 = length(il);         

%%-----  run opf  -----
gh_fcn = @(x)opf_consfcn(x, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);
hess_fcn = @(x, lambda, cost_mult)opf_hessfcn(x, lambda, cost_mult, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);

