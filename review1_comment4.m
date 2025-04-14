%% two_bus_opf_feasible_region.m
% This script implements the 2-bus OPF example and highlights the feasible
% region in the (V1, V2) plane where all constraints are simultaneously
% satisfied (P_G1, Q_G1, S12^2, S21^2, voltage bounds, etc.).
%
% The final plot shows:
%   - A gray-shaded area for the feasible set.
%   - Colored lines representing each constraint boundary.
%   - Axis ranges adjusted to reveal both feasible and infeasible areas.

clear; close all; clc;

%% 1) System Parameters

% Line impedance
z = 0.04 + 1i*0.2;   
abs_z = abs(z);   % |z|

% Admittance
y = 1/z;
G = real(y);
B = imag(y);

% Elements of the Ybus matrix (2-bus):
%   Y = [ y   -y
%        -y    y ]
G11 =  G;    B11 =  B;
G12 = -G;    B12 = -B;
G21 =  G12;  B21 =  B12;
G22 =  G11;  B22 =  B11;

%% 2) Load and Generator Data

% Active/Reactive load at bus 2
Pl2 = 3.5;          % P_load2
Ql2 = -3.5;          % Q_load2 (the paper uses Ql2 = -3.5 as injection sign)

% Generator limits at bus 1
PG_min = 0;   PG_max = 5;
QG_min = -4;  QG_max = 4;

% Squared line-flow limits
S_limit_sq = 30;   % for both S12^2 and S21^2

% Although the paper mentions Vmin=0.95, Vmax=1.05 for each bus,
% we will scan a *wider* range so we can see the boundary areas.
V1_min_plot = 0.60;   V1_max_plot = 1.20;
V2_min_plot = 0.60;   V2_max_plot = 1.20;

%% 3) Build Voltage Grid and Compute Power Flows

nPoints = 300;  % number of points for each voltage dimension
V1_vec = linspace(V1_min_plot, V1_max_plot, nPoints);
V2_vec = linspace(V2_min_plot, V2_max_plot, nPoints);

[V1, V2] = meshgrid(V1_vec, V2_vec);

%--- From bus-2 power balance, compute theta2 = atan2(N, D) ---

N = ( (-Pl2 - (V2.^2)*G22).*B21 ) - ( (-Ql2 + (V2.^2)*B22).*G21 );
D = ( (-Pl2 - (V2.^2)*G22).*G21 ) - ( (-Ql2 + (V2.^2)*B22).*B21 );
theta2 = atan2(N, D);

%--- Compute PG1, QG1 from bus-1 balance equations ---
PG1 = (V1.^2)*G11 + (V1.*V2).*( G12.*cos(theta2) - B12.*sin(theta2) );
QG1 = -(B11)*(V1.^2) + (V1.*V2).*( -G12.*sin(theta2) - B12.*cos(theta2) );

%--- Compute S12^2 and S21^2 ---
term = (V1.^2 + V2.^2 - 2.*V1.*V2.*cos(theta2));
S12_sq = ( (V1.^2) .* term ) / (abs_z^2);
S21_sq = ( (V2.^2) .* term ) / (abs_z^2);

%% 4) Construct the Feasible Region Mask

% We want the set of all (V1, V2) points where:
%   0.95 <= V1 <= 1.05
%   0.95 <= V2 <= 1.05
%   PG_min <= PG1 <= PG_max
%   QG_min <= QG1 <= QG_max
%   S12_sq <= S_limit_sq
%   S21_sq <= S_limit_sq

feasibleMask = ...
    (V1 >= 0.95) & (V1 <= 1.05) & ...
    (V2 >= 0.95) & (V2 <= 1.05) & ...
    (PG1 >= PG_min) & (PG1 <= PG_max) & ...
    (QG1 >= QG_min) & (QG1 <= QG_max) & ...
    (S12_sq <= S_limit_sq) & ...
    (S21_sq <= S_limit_sq);

% Convert logical mask to 0/1 for contour plotting
feasibleDouble = double(feasibleMask);

%% 5) Plot the Feasible Region + Constraint Boundaries

figure('Name','Feasible Region for 2-Bus OPF','Color','w');
hold on; box on;

% --- First, plot the "feasibility" as a filled contour:
% We contour at a level between 0 and 1 so that 1 => feasible region is shaded.
contourf(V1, V2, feasibleDouble, [0.5 1.5], ...
         'LineColor','none');
colormap([1 1 1; 0.6 0.6 0.6]);    % white for infeasible, gray for feasible
caxis([0 1]);                     % fix color scale

% --- Next, overlay boundary lines for each constraint. ---
% (You can change colors/styles to taste.)

% 1) PG1 = PG_min and PG1 = PG_max
contour(V1, V2, PG1, [PG_min PG_min], 'b--','LineWidth',1.5); 
contour(V1, V2, PG1, [PG_max PG_max], 'b-','LineWidth',1.5);

% 2) QG1 = QG_min and QG1 = QG_max
contour(V1, V2, QG1, [QG_min QG_min], 'r--','LineWidth',1.5);  
contour(V1, V2, QG1, [QG_max QG_max], 'r-','LineWidth',1.5);

% 3) S12^2 = S_limit_sq and S21^2 = S_limit_sq
contour(V1, V2, S12_sq, [S_limit_sq S_limit_sq], 'g--','LineWidth',1.5);
contour(V1, V2, S21_sq, [S_limit_sq S_limit_sq], 'g-','LineWidth',1.5);

% 4) Voltage boundaries: V1=0.95, V1=1.05, V2=0.95, V2=1.05
%    We can plot these as vertical/horizontal lines
plot([0.95 0.95],[V2_min_plot V2_max_plot],'k-','LineWidth',1);
plot([1.05 1.05],[V2_min_plot V2_max_plot],'k-','LineWidth',1);
plot([V1_min_plot V1_max_plot],[0.95 0.95],'k-','LineWidth',1);
plot([V1_min_plot V1_max_plot],[1.05 1.05],'k-','LineWidth',1);

% Axis settings
xlabel('V_1 (p.u.)');
ylabel('V_2 (p.u.)');
xlim([V1_min_plot, V1_max_plot]);
ylim([V2_min_plot, V2_max_plot]);
title('Feasible Region of the 2-Bus OPF Problem');

hold off;


%% check the point
V1=0.95;V2=1.025;
N =  (-Pl2 - (V2.^2)*G22).*B21  -  (-Ql2 + (V2.^2)*B22).*G21 ;
D =  (-Pl2 - (V2.^2)*G22).*G21  -  (-Ql2 + (V2.^2)*B22).*B21 ;
theta2 = atan2(N, D)
%--- Compute PG1, QG1 from bus-1 balance equations ---
PG1 = (V1.^2)*G11 + (V1.*V2).*( G12.*cos(theta2) - B12.*sin(theta2) )
QG1 = -(B11)*(V1.^2) + (V1.*V2).*( -G12.*sin(theta2) - B12.*cos(theta2) )
%--- Compute S12^2 and S21^2 ---
term = (V1.^2 + V2.^2 - 2.*V1.*V2.*cos(theta2));
S12_sq = ( (V1.^2) .* term ) / (abs_z^2)
S21_sq = ( (V2.^2) .* term ) / (abs_z^2)

Pmis = (V2.*V1).*( G21.*cos(theta2) + B21.*sin(theta2) )+(V2.^2)*G22+Pl2
Qmis = (V2.*V1).*( G21.*sin(theta2) - B21.*cos(theta2) )-(V2.^2)*B22+Ql2

