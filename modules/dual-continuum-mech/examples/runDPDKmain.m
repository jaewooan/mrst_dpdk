%% Example of a poromechanical dual-continuum simulation on a 2D grid with
%  stiff fractures
%
% The test problem is a uniaxial single-phase consolidation problem, such that the
% bottom boundary is fixed, and the left, right and top boundaries admit
% only vertical displacement. All of the boundaries are no flux boundaries,
% apart from the top boundary which is drained. 
%
% Within this example we demonstrate the following solvers:
%   * 'fully coupled'          : fully coupled solver
clear
clc
close all

%% Load required modules
mrstModule add dual-continuum-mech ad-core ad-mechanics dual-porosity ad-props vemmech

%% Input setting
nx = 10; ny = 10; 
Lx = 100; Ly = 100; %m
km = 1; kf = 1; %md
porom = 0.375; porof = 0.375;
vol_f = 0.0001; % frac vol fraction
mu = 0.0981; % centipoise = 1e-3 pa
rho = 1000; %kg/m3
cf = 4.4e-10; %fluid compressibility 1/pa
pref = 0; %pa
nu = 0.25; % Poisson's ratio
E_m = 1e9; % Young's modulus of matrix continuum
nu_m = nu; % Poisson's ratio of matrix continuum
E_f = 1e9; % Young's modulus of matrix continuum
nu_f = nu; % Poisson's ratio of fracture continuum
K_s = -999; % Solid stiffness, for biot = 1; K/(1-biot_coeff);%&70E9; % Solid stiffness
K_m = E_m / 3 / (1 - 2 * nu_m);
K_f = E_f / 3 / (1 - 2 * nu_f);
G_shear = 0.5 * E_m / (1 + nu_m);
if K_s < 0
    biot_f = 1;
    biot_m = 1;
    M_m = 1 / (cf*porom);
    M_f = 1 / (cf*porof);
else
    biot_f = 1 - K_f /K_s;
    biot_m = 1 - K_m /K_s;
    M_m = 1 / (cf*porom) + (biot_m-poro_m)/K_s;
    M_f = 1 / (cf*porof) + (biot_f-poro_f)/K_s;
end
B_m = biot_m * M_m / (K_m + biot_m^2*M_m);
B_f = biot_f * M_f / (K_f + biot_f^2*M_f);
nu_mu = (3*nu_m+biot_m* B_m*(1-2*nu_m))/(3-biot_m*B_m*(1-2*nu_m));
nu_fu = (3*nu_f+biot_f* B_f*(1-2*nu_f))/(3-biot_f*B_f*(1-2*nu_f));
force = 1e7; %pa boundary force
pres_m = 1/3*B_m*(1+nu_mu)*force;
pres_f = 1/3*B_f*(1+nu_fu)*force;


%% Setup default options
opt = struct('cartDims'            , [nx, ny], ...
             'L'                  , [Lx, Ly], ...
             'fluid_model'        , 'water', ...
             'verbose'            , false);

             
%% Setup Grid
G = cartGrid(opt.cartDims, opt.L);
G = computeGeometry(G);
%plotGrid(G);


%% Setup rock parameters (for flow)
% fracture
perm_fracture = kf*milli*darcy*ones(G.cartDims); 
rock_fracture = struct('perm', reshape(perm_fracture, [], 1), ...
              'poro', ones(G.cells.num, 1)*porom,... 
              'vol_fraction', ones(G.cells.num, 1)*vol_f); 
               % note, poro = intrinsic_poro*vol_fraction, therefore, for
               % this case intrinsic_poro is 0.6

% matrix
perm_matrix = km*milli*darcy*ones(G.cartDims); 
rock_matrix = struct('perm', reshape(perm_matrix, [], 1), ...
              'poro', ones(G.cells.num, 1)*porof,...
              'vol_fraction', ones(G.cells.num, 1)-rock_fracture.vol_fraction);
          
          
%% Setup fluid model, the current implementation only admits single-phase flow
% but the code is sufficiently generalisable to admit multi-phase flow. 
fluid = initSimpleADIFluid('phases', 'W', 'mu', mu*centi*poise, 'rho', ...
                                   rho*kilogram/meter^3, 'c', ...
                                   cf, 'pRef', pref);
fluid_fracture = fluid;
fluid_matrix = fluid; 


%% Setup material parameters for Biot and mechanics
E = HS_bound(E_m, E_f, nu_m, nu_f, rock_matrix.vol_fraction(1),...
             rock_fracture.vol_fraction(1), 'lower'); % Young's modulus of fractured rock mass with HS lower bound


E = repmat(E, G.cells.num, 1);
nu = repmat(nu, G.cells.num, 1);
E_m = repmat(E_m, G.cells.num, 1);
nu_m = repmat(nu_m, G.cells.num, 1);
E_f = repmat(E_f, G.cells.num, 1);
nu_f = repmat(nu_f, G.cells.num, 1);
K_s = repmat(K_s, G.cells.num, 1);


%% Setup boundary conditions for mechanics
% we first want to create a structure 'bc', which we can fudge by
% initialising the bc's using pside. 
oside = {'WEST', 'EAST', 'SOUTH', 'NORTH'};
bc = cell(4,1);
% for i = 1:numel(oside)
%     bc{i} = pside([], G, oside{i}, 0);
%     bc{i} = rmfield(bc{i}, 'type'); 
%     bc{i} = rmfield(bc{i}, 'sat');    
% end

for i = 1:numel(oside)
    if i == 2
        bc{i} = pside([], G, oside{i}, 0);
    else
        bc{i} = fluxside([], G, oside{i}, 0, 'sat', 1);
    end
    bc{i} = rmfield(bc{i}, 'type'); 
    bc{i} = rmfield(bc{i}, 'sat');    
end

% Displacement BCs
% Find the nodes for the different sides and set the boundaray conditions for
% elasticity.
for i = 1 : 4
    inodes = mcolon(G.faces.nodePos(bc{i}.face), G.faces.nodePos(bc{i}.face + 1) - 1);
    nodes = unique(G.faces.nodes(inodes));
    disp_bc = struct('nodes'   , nodes,      ...
                     'uu'      , 0,          ...
                     'faces'   , bc{i}.face, ...
                     'uu_face' , 0,          ...          
                     'mask'    , true(numel(nodes), G.griddim));
    bc{i}.el_bc = struct('disp_bc', disp_bc, 'force_bc', []);
end
bcdisp_zero = @(x) x*0.0; % Boundary displacement function set to zero.
bc_el_sides{1} = bc{1}; 
bc_el_sides{1}.el_bc.disp_bc.mask(:, 2) = false;   % west x fixed, y free
bc_el_sides{2} = bc{2}; 
bc_el_sides{2}.el_bc.disp_bc.mask(:, :) = false;   % east x free, y free
bc_el_sides{3} = bc{3}; 
bc_el_sides{3}.el_bc.disp_bc.mask(:, 1) = false;    % south x free, y fixed
bc_el_sides{4} = bc{4}; 
bc_el_sides{4}.el_bc.disp_bc.mask(:, :) = false;   % north x free, y fixed

% collect the displacement boundary conditions
nodes = [];
faces = [];
mask = [];
for i = 1 : numel(bc)
    if(~isempty(bc_el_sides{i}))
        nodes = [nodes; bc_el_sides{i}.el_bc.disp_bc.nodes]; %#ok
        faces = [faces; bc_el_sides{i}.el_bc.disp_bc.faces]; %#ok
        mask  = [mask; bc_el_sides{i}.el_bc.disp_bc.mask]; %#ok
    end
end

disp_node = bcdisp_zero(G.nodes.coords(nodes, :));
disp_faces = bcdisp_zero(G.faces.centroids(faces, :)); 
disp_bc = struct('nodes', nodes, 'uu', disp_node, 'faces', faces, 'uu_face', disp_faces, 'mask', mask); 

% Force BCs
facesf = bc{4}.face;
force_bc = struct('faces', facesf, 'force', force*ones(G.cartDims(1),1)*[0 -1]);
el_bc = struct('disp_bc', disp_bc, 'force_bc', force_bc);


%% Gravity
% the gravity in this option affects only the fluid behaviour
gravity off;
    

%% Setup load for mechanics
% in this example we do not impose any volumetric force
load = @(x) (0*x);


%% Gather all the mechanical parameters in a struct that will be used to
% to define the mech problem
mech = struct('E', E, 'nu', nu, 'E_m', E_m, 'nu_m', nu_m, 'E_f', E_f, 'nu_f', nu_f,...
              'K_s', K_s, 'el_bc', el_bc, 'load', load);

%% Setup fully coupled and fixed stress splitting models
fullycoupledOptions = {'verbose', opt.verbose};
DC_model = DualContMechWaterModel(G, {rock_fracture, rock_matrix}, {fluid_fracture, fluid_matrix}, mech, fullycoupledOptions{:});
DC_model.nonlinearTolerance = 1e-8;
d1 = 0.1; % spacing of fracture set 1
d2 = d1; % spacing of fracture set 2
fracture_spacing = repmat([d1,d2],G.cells.num,1);
shape_factor_name = 'Lim_AzizShapeFactor';
DC_model.transfer_model_object = SimpleTransferFunction(shape_factor_name, fracture_spacing);
DC_model = DC_model.validateModel();

%% Setup initial state and fluid BCs
pressure = ones(G.cells.num,1)*44.17e5; % pa
u_init = zeros(2*(nx+1)*(ny+1), 1);
u_init(1:2:end) = force*nu_fu/2/G_shear*G.nodes.coords(:,1);
u_init(2:2:end) = -force*(1-nu_fu)/2/G_shear*G.nodes.coords(:,2);

% state0 = struct('pressure', pressure, 'pressure_matrix', pressure, 's', ones(G.cells.num, 1), 'swm', ones(G.cells.num, 1));
state0 = struct('pressure', ones(G.cells.num,1)*pres_f, 'pressure_matrix', ones(G.cells.num,1)*pres_m, 's', ones(G.cells.num, 1), 'swm', ones(G.cells.num, 1),...
                'u', u_init);

% need to initiate the fluid bc's, bc's are the same for micro and macro scales 
%bc_f0 = pside([], G, 'WEST', 10000000, 'sat', 1);
bc_f0 = fluxside([], G, 'WEST', 0, 'sat', 1);
%bc_f0 = fluxside(bc_f0, G, 'EAST', 0,'sat', 1);
bc_f0 = pside(bc_f0 , G, 'EAST', 0, 'sat', 1);
bc_f0 = fluxside(bc_f0, G, 'SOUTH', 0, 'sat', 1);
bc_f0 = fluxside(bc_f0, G, 'NORTH', 0, 'sat', 1);
%bc_f0 = pside(bc_f0 , G, 'SOUTH', 0, 'sat', 1);
%bc_f0 = pside(bc_f0 , G, 'NORTH', 0, 'sat', 1);


%% Simulate 
time = [0,linspace(0.1148543e-02,9*0.1148543e-02, 9),linspace(0.1148543e-01,9*0.1148543e-01, 9),...
       linspace(0.1148543,9*0.1148543, 9), linspace(0.1148543e1,9*0.1148543e1, 9),...
       linspace(0.1148543e2,10*0.1148543e2, 10)]*3600*24;
% dt = 1;
% time = [0, linspace(dt, 10, 10*dt)];
dt = diff(time);
[p_m, p_f, u, states] = simDC_mech(state0, dt, DC_model, bc_f0);


%% Plot results
figure
semilogx(time, p_m(1, :), '-', 'linewidth', 1.5)
hold on
semilogx(time, p_f(1, :), '--', 'linewidth', 1.5, 'markersize', 6)
hold on
xlabel('time [s]')
ylabel('average pressure [Pa]')
legend('matrix', 'fracture')
title('Results for the intrinsic fracture stiffness simulation')

%% test
close all
figure
for i=1:nx
    X_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,1)';
    Y_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,2)';
    R_(i,:) = p_f(nx*(i-1)+1:nx*i,size(time,2));
end
pcolor(X_, Y_, R_);
shading interp
set(gca,'YDir','normal') 
colorbar

figure
for i=1:nx
    X_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,1)';
    Y_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,2)';
    R_(i,:) = p_m(nx*(i-1)+1:nx*i,size(time,2));
end
pcolor(X_, Y_, R_);
shading interp
set(gca,'YDir','normal') 
colorbar


figure
for i=1:nx
    X_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,1)';
    Y_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,2)';
    R_(i, :) = states{1}.stress(nx*(i-1)+1:nx*i, 2);
end
pcolor(X_, Y_, R_);
shading interp
set(gca,'YDir','normal') 
colorbar

%%
figure
for i=1:nx+1
    X__(i, :) =  G.nodes.coords((nx+1)*(i-1)+1:(nx+1)*i,1)';
    Y__(i, :) =  G.nodes.coords((nx+1)*(i-1)+1:(nx+1)*i,2)';
    R__(i,:) = u((nx+1)*2*(i-1)+1:2:(nx+1)*2*i,size(time,2));
end
pcolor(X__, Y__, R__);
set(gca,'YDir','normal') 
colorbar


figure
for i=1:nx+1
    X__(i, :) =  G.nodes.coords((nx+1)*(i-1)+1:(nx+1)*i,1)';
    Y__(i, :) =  G.nodes.coords((nx+1)*(i-1)+1:(nx+1)*i,2)';
    R__(i, :) = u((nx+1)*2*(i-1)+2:2:(nx+1)*2*i+1,5);
end
pcolor(X__, Y__, R__);
set(gca,'YDir','normal') 
colorbar



%% Mandel - analytical solution 
biot_coefficient = 1;
biot_modulus = 1/(porom*cf); % M

sol = MandelAnalyticalSolution(Lx, Ly, E(1), nu(1), biot_coefficient(1), biot_modulus(1), km*milli*darcy, mu*centi*poise, force);
dx = Lx/nx;
dy = Ly/ny;
pcoord = [dx/2; dy/2];
ucoord = [0; Ly];

% compute analytical solution for pressure at x=dx, z = dx (center of first cell),
%                          for displacement at x =0, z = 100

for i = 1:length(time)
  t = time(i);
  sol = compute_pressure_solution(sol, t, pcoord);
  sol = compute_displacement_solution(sol, t, ucoord);
  p(i) = sol.pressure;
  u_anal(i) = abs(sol.displacement(2));
end
maxt = max(time);

% pressure at first cell (bottom left) 
figure; set(gcf,'color','white'); set(gcf, 'Position', [400 200 800 600])
% semilogx(time, p/p(1), 'r', ...
%                  time, p_m(1, :)/p_m(1, 1), 'o', ...
%                  time, p_f(1, :)/p_f(1, 1), 'k*');
semilogx(time, p, 'r', ...
                 time, p_m(1, :), 'o', ...
                 time, p_f(1, :), 'k*');
xlabel('Time, Days'); ylabel('Normalized Pressure');
leg = legend('analytical', 'Matrix', 'Fracture', 'Location', 'SouthWest');
pos = get(leg, 'Position');
pos = [pos(1), pos(2)*1.5, pos(3), pos(4)];
set(leg, 'Position', pos);
% ylim([0, 1.2]);     xlim([0, maxt]);
set(gca,'XTick',[1e-2 1e-1 1 1e1]);
%Image = getframe(gcf); %imwrite(Image.cdata, 'mandel_pressure.png');

% displacement at x =0, z = 102
figure; set(gcf,'color','white'); set(gcf, 'Position', [400 200 800 600])
% semilogx(time, u_anal/u_anal(1), 'r', time, u(421*2,:)/u(421*2,1), 'o');
semilogx(time, -u_anal, 'r', time, u(end,:), 'o');
xlabel('Time, Days'); ylabel('Normalized Settlement');
leg = legend('analytical', 'y Displacment', 'Location', 'SouthWest');
pos = get(leg, 'Position');
pos = [pos(1), pos(2)*1.5, pos(3), pos(4)];
set(leg, 'Position', pos);
% xlim([0, maxt]); ylim([0.99, 1.5]);
set(gca,'XTick',[1e-2 1e-1 1 1e1]);
%Image = getframe(gcf); imwrite(Image.cdata, 'mandel_displacement.png');
 


        