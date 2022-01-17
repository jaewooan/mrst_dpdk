%% Workflow example for the MICP model in a 3D system assessing CO2 leakage
%
% This example aims to show complete workflow for creating, running, and
% analyzing a 3D-flow system with a leakage path using the MICP
% mathematical model. To asses the CO2 distribution before and after
% treatment, we use the 'TwoPhaseWaterGasModel' in the MRST co2lab module
% (we use the system relationships/properties as in the 'basic_3D_example'
% in the co2lab module).
%
% For details on the MICP model, see
% Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021. Practical
% approaches to study microbially induced calcite precipitation at the
% field scale. Int. J. Greenh. Gas Control 106, 103256.
% https://doi.org/10.1016/j.ijggc.2021.103256

% To get distmesh for first time, uncomment and run the 4 following lines
% pth = fullfile(ROOTDIR, 'utils', '3rdparty', 'distmesh');
% mkdir(pth)
% unzip('http://persson.berkeley.edu/distmesh/distmesh.zip', pth);
% mrstPath('reregister', 'distmesh', pth);

% Required modules
pth = fullfile(ROOTDIR, 'utils', '3rdparty', 'distmesh');
mrstPath('reregister', 'distmesh', pth);
mrstModule add ad-blackoil ad-core ad-micp ad-props co2lab distmesh

%% Reservoir geometry/properties and model parameters
%
% The domain has a length of 200 m, a height of 160 m, and width of 180 m.
% We remove the domain cells where there is caprock. The leakage
% xy-position is set to the origin (0, 0). The grid is coarser in order to
% run the example in couple of minutes. You should try with larger domains
% and finer grids (see https://doi.org/10.1016/j.ijggc.2021.103256).

L = 200;            % Reservoir length, m
W = 180;            % Reservoir width, m
ht = 30;            % Top aquifer heigth, m
hl = 100;           % Leakage heigth, m
hb = 30;            % Bottom aquifer heigth, m
H = ht + hl + hb;   % Reservoir heigth, m
D = 1500;           % Depth of aquifer top surface, m
a = 1;              % Leakage width, m
w = -20;            % Well x-location from leakage, m
hmin = a;           % Minimum grid size, m
hmid = 3;           % Medium grid size, m
hmax = L;           % Maximum grid size, m
B = 10;             % hmid to hmax transition radius, m
% Here we create the mesh using distmesh. For information/tutorials, visit
% http://persson.berkeley.edu/distmesh/
fd = @(p) drectangle(p, -L / 2, L / 2, -W / 2, W / 2);
fh = @(p) min(min(hmid + 0.3 * abs(dcircle(p, w, 0, 0)), hmid) .* ...
 (abs(dcircle(p, w, 0, 0)) < B) + min(hmid + 0.3 * abs(dcircle(p, w, 0, ...
        B)), hmax) .* (abs(dcircle(p, w, 0, 0)) >= B), min(hmin + 0.3 * ...
         abs(dcircle(p, 0, 0, 0)), hmid) .* (abs(dcircle(p, 0, 0, 0)) < ...
                  B) + min(hmid + 0.3 * abs(dcircle(p, 0, 0, B)), hmax) ...
                                       .* (abs(dcircle(p, 0, 0, 0)) >= B));
[p, t] = distmesh2d(fd, fh, hmin, [-L / 2, -W / 2; L / 2, W / 2], ...
       [-L / 2, -W / 2; L / 2, -W / 2; -L / 2, W / 2; L / 2, W / 2; 0, 0]);
close
Z = [0 3 10 ht : 10 : ht + hl ht + hl + 3 ht + hl + 10 H];
G = makeLayeredGrid(pebi(triangleGrid(p, t)), max(size(Z)) - 1);
zn = G.nodes.num / max(size(Z));
for i = 1 : max(size(Z))
    G.nodes.coords(1 + zn * (i - 1) : 1 : zn * i, 3) = Z(i) + D;
end
G = computeGeometry(G);
c = G.cells.centroids;
G = removeCells(G, c(:, 1) .^ 2 + c(:, 2) .^ 2 > (hmin / 2) ^ 2 & ...
                                 c(:, 3) < D + ht + hl & c(:, 3) > D + ht);
G = computeGeometry(G);
c = G.cells.centroids;
C = ones(G.cells.num, 1);

% Rock
K0 = 2e-14 * C;              % Aquifer permeability, m^2
% We find the cells where the leakage path is located and set a higher
% permeability.
cellsl = G.cells.indexMap;
cellsl = cellsl(c(:, 1) .^ 2 + c(:, 2) .^ 2 < (hmin / 2) ^ 2 & ...
                                  c(:, 3) < D + H - hb & c(:, 3) > D + ht);
cellsF =  G.cells.indexMap;
idx = ismember(cellsF, cellsl);
K0(idx) = 1e-12;             % Leakage permeability, m^2
porosity = 0.15;             % Aquifer/leakage porosity, [-]
rock = makeRock(G, K0, porosity);

% Fluid properties
fluid.muw = 2.535e-4;        % Water viscocity, Pa s
fluid.bW = @(p) 0 * p + 1;   % Water formation volume factor, [-]
fluid.rhoWS = 1045;          % Water density, kg/m^3

% Remaining model parameters (we put them on the fluid structure)
fluid.rho_b = 35;            % Density (biofilm), kg/m^3
fluid.rho_c = 2710;          % Density (calcite), kg/m^3
fluid.k_str = 2.6e-10;       % Detachment rate, m/(Pa s)
fluid.diffm = 2.1e-9;        % Diffusion coefficient (microbes), m^2/s
fluid.diffo = 2.32e-9;       % Diffusion coefficient (oxygen), m^2/s
fluid.diffu = 1.38e-9;       % Diffusion coefficient (urea), m^2/s
fluid.alphaL = 1e-3;         % Disperison coefficient (longitudinal), m
fluid.alphaT = 4e-4;         % Disperison coefficient (transverse), m
fluid.eta = 3;               % Fitting factor, [-]
fluid.k_o = 2e-5;            % Half-velocity constant (oxygen), kg/m^3
fluid.k_u = 21.3;            % Half-velocity constant (urea), kg/m^3
fluid.mu = 4.17e-5;          % Maximum specific growth rate, 1/s
fluid.mu_u = 0.0161;         % Maximum rate of urease utilization, 1/s
fluid.k_a = 8.51e-7;         % Microbial attachment rate, 1/s
fluid.k_d = 3.18e-7;         % Microbial death rate, 1/s
fluid.Y = 0.5;               % Yield growth coefficient, [-]
fluid.Yuc = 1.67;            % Yield coeccifient (calcite/urea), [-]
fluid.F = 0.5;               % Oxygen consumption factor, [-]
fluid.crit = 0.1;            % Critical porosity, [-]
fluid.kmin = 1e-20;          % Minimum permeability, m^2
fluid.cells = C;             % Array with all cells, [-]
fluid.ptol = 1e-4;           % Porosity tolerance to stop the simulation

% Porosity-permeability relationship
fluid.K = @(poro) (K0 .* ((poro - fluid.crit) / (porosity - fluid.crit))...
               .^ fluid.eta + fluid.kmin) .* K0 ./ (K0 + fluid.kmin) .* ...
                  (poro > fluid.crit) + fluid.kmin .* (poro <= fluid.crit);

% The two following lines are not really used in these simulations since
% the current MICP implementation only considers single-phase flow (it is
% possible to extend to two-phase flow), but since the implementation is
% based on the 'equationsOilWaterPolymer' script (two-phase flow), they are
% required to avoid errors.
fluid.bO   = fluid.bW;
fluid.rhoOS = fluid.rhoWS;

% Gravity
gravity on

%% CO2 assesment
%
% We simulate the CO2 distribution on the domain before MICP treatment.
% We use the 'TwoPhaseWaterGasModel' in the MRST co2lab module (the system
% relationships/properties are set as in the 'basic_3D_example' in the
% co2lab module.

state0.pressure = fluid.rhoWS * norm(gravity) * c(:, 3); % Initial pressure
state0.s = repmat([1, 0], G.cells.num, 1); % Initial saturations
co2     = CO2props(); % Load sampled tables of co2 fluid properties
p_ref   = 15 * mega * Pascal; % Reference pressure
t_ref   = 70 + 273.15; % Reference temperature, in Kelvin
rhoco2  = co2.rho(p_ref, t_ref); % CO2 density at ref. press/temp
cf_co2  = co2.rhoDP(p_ref, t_ref) / rhoco2; % CO2 compressibility
cf_wat  = 0; % Water compressibility (zero)
cf_rock = 0; % Rock compressibility (zero)
muco2   = co2.mu(p_ref, t_ref) * Pascal * second; % CO2 viscosity

% Use function 'initSimpleADIFluid' to make a simple fluid object
fluidH2OCO2 = initSimpleADIFluid('phases', 'WG'           , ...
                           'mu'  , [fluid.muw, muco2]     , ...
                           'rho' , [fluid.rhoWS, rhoco2]  , ...
                           'pRef', p_ref                  , ...
                           'c'   , [cf_wat, cf_co2]       , ...
                           'cR'  , cf_rock                , ...
                           'n'   , [2 2]);

% Change relperm curves
srw = 0.27;
src = 0;
fluidH2OCO2.krW = @(s) fluidH2OCO2.krW(max((s - srw) ./ (1 - srw), 0));
fluidH2OCO2.krG = @(s) fluidH2OCO2.krG(max((s - src) ./ (1 - src), 0));

% Add capillary pressure curve
pe = 5 * kilo * Pascal;
pcWG = @(sw) pe * sw .^ (-1 / 2);
fluidH2OCO2.pcWG = @(sg) pcWG(max((1 - sg - srw) ./ (1 - srw), 1e-5));

% Create Model
model = TwoPhaseWaterGasModel(G, rock, fluidH2OCO2, 0, 0);

% Create Well
QCO2 = 80 / day;      % Injection rate, m^3/day
r = 0.15;             % Well radius, m
[~, iw] = min(abs((c(:, 1) - w) .^ 2 + c(:, 2) .^ 2));
cellsW =  1 : G.cells.num;
cellsW = cellsW(abs(c(:, 1) - c(iw, 1)) < eps & ...
                     abs(c(:, 2) - c(iw, 2)) < eps & c(:, 3) > D + H - hb);
% Injector
WCO2 = addWell([], G, rock, cellsW, 'Type', 'rate', 'Comp_i', [0, 1], ...
                                                 'Val', QCO2, 'Radius', r);

% We set a constant hydrostatic pressure on the left and right side of the
% domain to model open boundaries.
f = boundaryFaces(G);
f = f(abs(G.faces.normals(f, 1)) > eps & (G.faces.centroids(f, 1) < ...
                  -hmax / 4 | G.faces.centroids(f, 1) > L / 2 - hmax / 4));
fp = G.faces.centroids(f, 3) * fluid.rhoWS * norm(gravity);
bc = addBC([], f, 'pressure', fp, 'sat', [0 0]);

% Setup some schedule
dt_co2 = day;
nt = 30 * day / dt_co2;
timestepsCO2 = repmat(dt_co2, nt, 1);

% Make schedule
schedule = simpleSchedule(timestepsCO2, 'W', WCO2, 'bc', bc);

% Simulate
if exist('OCTAVE_VERSION', 'builtin') == 0
  % The two last entries in the 'getPlotAfterStepCO2' function are the
  % azimuth and elevation angles for view of the current axes while
  % visualizing the solution. This function could be also used in
  % other mrst examples to visualize the evolution of the non-wetting phase
  % saturation, e.g., it works in the 'example3D' in the co2lab module, if
  % line 145 is replaced by the three following lines:
  % mrstModule add ad-micp
  % fn = getPlotAfterStepCO2(initState, model, 340, 20);
  % [wellSol, states] = simulateScheduleAD(initState, model, schedule, ...
  % 'afterStepFn', fn);
  fn = getPlotAfterStepCO2(state0, model, 340, 20);
  [~, statesCO2beforeMICP] = simulateScheduleAD(state0, model, ...
                                              schedule, 'afterStepFn', fn);
else
  [~, statesCO2beforeMICP] = simulateScheduleAD(state0, model, schedule);
end

% Compute leakage rate (CO2 rate through the lowest grid face on the
% leakage path). If the grid properties are modified, then it could be
% necessary to also modify the choosen value 0.01 in L 225-226.
lrbeforeMICP = zeros(nt, 1);
fc = G.faces.centroids;
facel =  1 : G.faces.num;
facel = facel(fc(:, 3) < D + hl + ht + 0.01 & fc(:, 3) > D + ...
        hl + ht - 0.01 & (fc(:, 1) .^ 2 + fc(:, 2) .^ 2 < (hmin / 2) ^ 2));
for i = 1 : nt
    lrbeforeMICP(i) = abs(statesCO2beforeMICP{i}.flux(facel, 2));
end

%% MICP treatment

% Set the injection strategy.
N = 9; % Number of injection phases in the injection strategy
M = zeros(N, 6); % The entries per row are: time, dt, rate, m, o, and u.
dt_on = hour; % Time step when the well is on
dt_off = 5 * hour;  % Time step when the well is off

M(1, :) = [15 * hour,   dt_on,  0.3, 0.01, 0,      0];
M(2, :) = [11 * hour,   dt_on,  0.3, 0,    0,      0];
M(3, :) = [75 * hour,   dt_off, 0,   0,    0,      0];
M(4, :) = [30 * hour,   dt_on,  0.3, 0,    0.04,   0];
M(5, :) = [5 * hour,    dt_on,  0.3, 0,    0,      0];
M(6, :) = [25 * hour,   dt_off, 0,   0,    0,      0];
M(7, :) = [40 * hour,   dt_on,  0.3, 0,    0,     60];
M(8, :) = [10 * hour,   dt_on,  0.3, 0,    0,      0];
M(9, :) = [40 * hour,   dt_off, 0,   0,    0,      0];

% Create Well. For the injection of the MICP components, we create a well
% with an upper and lower parts, where the components are injected in the
% top part and only water on the lower part (hb is the bottom aquifer
% heigth).
Whu = 3 / hb;   % Upper fraction part of the well to inject the components
Whb = 1 - Whu;  % Bottom fraction part of the well to inject the components
cellsW =  1 : G.cells.num;
cellsWu = cellsW(abs(c(:, 1) - c(iw, 1)) < eps & abs(c(:, 2)-c(iw, 2)) ...
                < eps & c(:, 3) > D + H - hb & c(:, 3) < D + H - Whb * hb);
% Upper injector
W = addWell([],G, rock, cellsWu, 'Type', 'rate', 'Comp_i', [1, 0], ...
                                        'Val', Whu * M(1, 3), 'Radius', r);
cellsWb = cellsW(abs(c(:, 1) - c(iw, 1)) < eps & ...
               abs(c(:, 2) - c(iw, 2)) < eps & c(:, 3) > D + H - Whb * hb);
% Lower injector
W = addWell(W, G, rock, cellsWb, 'Type', 'rate', 'Comp_i', [1, 0], ...
                                        'Val', Whb * M(1, 3), 'Radius', r);

% Add the fields to the wells/bc for the additional components
bc.m = zeros(size(bc.sat, 1), 1);
bc.o = zeros(size(bc.sat, 1), 1);
bc.u = zeros(size(bc.sat, 1), 1);
bc.b = zeros(size(bc.sat, 1), 1);
bc.c = zeros(size(bc.sat, 1), 1);
for i = 1 : 2
    W(i).o = 0;
    W(i).u = 0;
    W(i).m = 0;
end
W(1).m = M(1, 4);  % Injected microbial concentration, kg/m^3
W(1).o = M(1, 5);  % Injected oxygen concentration, kg/m^3
W(1).u = M(1, 6);  % Injected urea concentration, kg/m^3
% The injection well is not on the boundary
G.injectionwellonboundary = 0;

% Setup some schedule
nt = sum(M(:, 1) ./ M(:, 2));
timesteps = repmat(dt_on, nt, 1);
schedule = simpleSchedule(timesteps, 'W', W, 'bc', bc);
for i = 2 : N
    schedule.control(i) = schedule.control(i - 1);
    schedule.step.control(sum(M(1 : i - 1, 1) ./ M(1 : i - 1, 2)) + 1 : ...
                                                                  end) = i;
    schedule.step.val(sum(M(1 : i - 1, 1) ./ M(1 : i - 1, 2)) + 1 : ...
                                                            end) = M(i, 2);
    schedule.control(i).W(1).val = Whu * M(i, 3);
    schedule.control(i).W(2).val = Whb * M(i, 3);
    schedule.control(i).W(1).m = M(i, 4);
    schedule.control(i).W(1).o = M(i, 5);
    schedule.control(i).W(1).u = M(i, 6);
end

% Maximum injected oxygen and urea concentrations.
fluid.Comax = max(M(:, 5));
fluid.Cumax = max(M(:, 6));

% Create model
model = MICPModel(G, rock, fluid);

% Initial Condition
state0_micp = state0;
state0_micp.m = zeros(G.cells.num, 1);
state0_micp.o = zeros(G.cells.num, 1);
state0_micp.u = zeros(G.cells.num, 1);
state0_micp.b = zeros(G.cells.num, 1);
state0_micp.c = zeros(G.cells.num, 1);

% If MATLAB is used, we use the getPlotAfterStepMICP function to visualize
% the results at each time step.
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    ok = 'true';
    fn = checkCloggingMICP(ok);
else
    fn = getPlotAfterStepMICP(state0_micp, model, 340, 20);
end
[~,statesMICP] = simulateScheduleAD(state0_micp, model, schedule, ...
                                                        'afterStepFn', fn);

%% CO2 assesment after MICP treatment
%
% We simulate the CO2 distribution on the domain after MICP treatment.

% Compute porosity and permeability after MICP treatment
porosityafterMICP = porosity - statesMICP{end}.c - statesMICP{end}.b;
KafterMICP = fluid.K(porosityafterMICP);
rock = makeRock(G, KafterMICP, porosityafterMICP);

% Create model
model = TwoPhaseWaterGasModel(G, rock, fluidH2OCO2, 0, 0);

% Make schedule
schedule = simpleSchedule(timestepsCO2, 'W', WCO2, 'bc', bc);

% Simulate
if exist('OCTAVE_VERSION', 'builtin') == 0
  fn = getPlotAfterStepCO2(state0, model, 340, 20);
  [~, statesCO2afterMICP] = simulateScheduleAD(state0, model, schedule, ...
                                                        'afterStepFn', fn);
else
  [~, statesCO2afterMICP] = simulateScheduleAD(state0, model, schedule);
end

% Compute the CO2 leakage rate after MICP treatment
lrafterMICP = zeros(max(size(timestepsCO2)), 1);
for i = 1 : size(timestepsCO2)
    lrafterMICP(i) = abs(statesCO2afterMICP{i}.flux(facel, 2));
end

%% Process the data
%
% Plot the comparison of CO2 leakage before and after micp treatment. Using
% the default example setting, we observe the leakage is reduced after MICP
% treatment but still there is significant leakage. Then additional MICP
% treatments could be applied to reduce the leakage as reported in
% https://doi.org/10.1016/j.ijggc.2021.103256

figure;
hold on
plot((1 : size(timestepsCO2)) * dt_co2 / day, lrbeforeMICP * 100 /...
             QCO2, 'color', [1 0.2 0.2], 'LineWidth', 9, 'LineStyle', '-');
plot((1 : size(timestepsCO2)) * dt_co2 / day, lrafterMICP * 100 / ...
               QCO2, 'color', [1 0.5 0], 'LineWidth', 9, 'LineStyle', '-');
hold off
legend('Before MICP', 'After MICP', 'Location', 'southeast');
xlabel('Time [d]');
ylabel('CO2 leakage rate/injection rate [%]');
grid on

% If Octave is used, then the results are printed in vtk format to be
% visualized in Paraview and the 'return' command is executed as currently
% it is not possible to run 'plotToolbar' in Octave.

% Write the results to be read in ParaView (GNU Octave)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mkdir vtk_3DCase;
    cd vtk_3DCase;
    mrsttovtk(G, statesCO2afterMICP, 'statesCO2afterMICP', '%f');
    mrsttovtk(G, statesMICP, 'statesMICP', '%f');
    mrsttovtk(G, statesCO2beforeMICP, 'statesCO2beforeMICP', '%f');
    return
end

% If MATLAB is used, then the plotToolbar is used to show the results. For
% this is necessary to add the mrst-gui module.
mrstModule add mrst-gui
figure;
plotToolbar(G, statesMICP, 'field', 's:1', 'lockCaxis', true);
view(340, 20); axis tight; colorbar; caxis([0 1]);

%% Copyright notice
%{
Copyright 2021, NORCE Norwegian Research Centre AS, Computational
Geosciences and Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}
