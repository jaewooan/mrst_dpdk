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

%% Basic Simulation input setting
nx = 10; ny = 10; 
Lx = 100; Ly = 100; %m
porom = 0.3; porof = 0.5;
rho = 1000; %kg/m3
cf = 1e-10; %fluid compressibility 1/pa
pref = 0; %pa
K_s = -999; % Solid stiffness, for biot = 1, Ks = infty = K / (1-biot_coeff)
time = linspace(0.01, 1, 100)*day; % sec
force = 5e7; % pa, boundary traction
p_init = 3e7; % pa, initial pressure
sat_init = 1;
vol_f = 0.02; % frac vol fraction
PV = Lx*Ly *((1-vol_f)*porom + vol_f*porof);
q_inj = 0.5*PV/time(end)/10; % m/sec
q_prod = q_inj; % m/sec
well_radius = 0.1; %m
d1 = 0.03;%Lx/nx/100; % m, spacing of fracture set 1
d2 = d1; % m, spacing of fracture set 2

%% Uncertainty generation
% Expected uncertainty for future..
nu = 0.25; % Poisson's ratio
nu_m = nu*ones(nx*ny,1); % Poisson's ratio of matrix continuum
nu_f = nu*ones(nx*ny,1); % Poisson's ratio of fracture continuum
E_m = 1e9*ones(nx*ny,1); % Young's modulus of matrix continuum
E_f = 1e6*ones(nx*ny,1); % Young's modulus of fracture continuum
mu = 1; % centipoise = 1e-3 pa

% our target for uncertainty: perm
% km = 1; kf = 100; %md
perm_matrix = 0.01*milli*darcy*ones(nx,ny); 
perm_fracture = 1000*milli*darcy*ones(nx,ny); 

% Load SGS perm files
load('perm.mat');
results=zeros(0,0);
for i=1:nCase
    %% Setup default options
%     E_m = reshape(EM_tot(:,i), nx*ny, 1);
%     E_f = reshape(EF_tot(:,i), nx*ny, 1);
%     perm_matrix = reshape(perm_m_tot(:,i), nx, ny)*milli*darcy*0.01;
%     perm_fracture = reshape(perm_f_tot(:,i), nx, ny)*milli*darcy;
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
    rock_fracture = struct('perm', reshape(perm_fracture, [], 1), ...
                  'poro', ones(G.cells.num, 1)*porom,... 
                  'vol_fraction', ones(G.cells.num, 1)*vol_f); 
                   % note, poro = intrinsic_poro*vol_fraction, therefore, for
                   % this case intrinsic_poro is 0.6
    
    % matrix
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
    [E, dd1, dd2, nu] = HS_bound(E_m, E_f, nu_m, nu_f, rock_matrix.vol_fraction,...
                 rock_fracture.vol_fraction, 'lower'); % Young's modulus of fractured rock mass with HS lower bound
    
    
    %E = repmat(E, G.cells.num, 1);
    
    
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
        bc{i} = fluxside([], G, oside{i}, 0);
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
    bc_el_sides{1}.el_bc.disp_bc.mask(:, 2) = false;   % west: x fixed, y free
    bc_el_sides{2} = bc{2}; 
    bc_el_sides{2}.el_bc.disp_bc.mask(:, :) = false;   % east x free, y free
    bc_el_sides{3} = bc{3}; 
    bc_el_sides{3}.el_bc.disp_bc.mask(:, 1) = false;    % south x free, y fixed
    bc_el_sides{4} = bc{4}; 
    bc_el_sides{4}.el_bc.disp_bc.mask(:, :) = false;   % north x free, y free
    
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
    nEast = 2;
    nNorth = 4;
    facesf = [bc{nEast}.face; bc{nNorth}.face]; % same traction at the top and right
    force_vec = [force*ones(size(bc{nEast}.face, 1),1)*[-1 0];...
                 force*ones(size(bc{nNorth}.face, 1),1)*[0 -1]];
    force_bc = struct('faces', facesf, 'force', force_vec);
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
    DC_model.nonlinearTolerance = 1e-10;
    fracture_spacing = repmat([d1,d2],G.cells.num,1);
    shape_factor_name = 'Lim_AzizShapeFactor';
    DC_model.transfer_model_object = SimpleTransferFunction(shape_factor_name, fracture_spacing);
    DC_model = DC_model.validateModel();
    
    %% Setup initial state and fluid BCs
    pressure = ones(G.cells.num,1)*p_init; % pa
    state0 = struct('pressure', pressure, 'pressure_matrix', pressure, 's', sat_init*ones(G.cells.num, 1), 'swm', sat_init*ones(G.cells.num, 1));
    
    % need to initiate the fluid bc's, bc's are the same for micro and macro scales 
    bc_f0 = fluxside([], G, 'WEST', 0, 'sat', 1); % no flow boundary
    bc_f0 = fluxside(bc_f0 , G, 'EAST', 0, 'sat', 1); % no flow boundary
    bc_f0 = fluxside(bc_f0, G, 'SOUTH', 0, 'sat', 1); % no flow boundary
    bc_f0 = fluxside(bc_f0, G, 'NORTH', 0, 'sat', 1); % no flow boundary
    %bc_f0 = pside(bc_f0 , G, 'SOUTH', 0, 'sat', 1); % free boundary
    
    %% Simulate 
    dt = diff(time);
    cellCenters = [nx*(ny/2-1) + nx/2+3];
    W = [];
%     W = addWell(W, DC_model.G, DC_model.rock, 1, 'Type', 'rate', ...
%          'comp_i', [1, 0, 0], 'Name', ['INJ'], 'Val', q_inj, 'sign', 1, 'Radius', well_radius, 'Dir','z');
%     W = addWell(W, DC_model.G, DC_model.rock, cellCenters, 'Type', 'rate', ...
%          'Name', ['PROD'], 'Val', -q_prod, 'sign', -1, 'Radius', well_radius, 'Dir','x');
    W = addWell(W, DC_model.G, DC_model.rock, cellCenters, 'Type', 'bhp', ...
         'Name', ['PROD'], 'Val', p_init, 'sign', -1, 'Radius', well_radius, 'Dir','z');
    [p_m, p_f, u, states] = simDC_mech(state0, dt, DC_model, bc_f0, W);
    result_ = [p_m;p_f;u(1:2:2*(nx+1)*(ny+1),:);u(2:2:2*(nx+1)*(ny+1),:)]; %pm, pf, ux, uy
    results = [results,reshape(result_, size(result_,1)*size(result_,2),1)];
end

name = sprintf('pm_pf_ux_uy_nt%d_ncase_%d.mat', length(time), nCase);
save(name, 'results');

%% Plot results
close all
figure
semilogx(time, p_m(nx*(ny/2-1) + nx/2, :), '-', 'linewidth', 1.5)
hold on
semilogx(time, p_f(nx*(ny/2-1) + nx/2, :), '--', 'linewidth', 1.5, 'markersize', 6)
hold on
xlabel('time [s]')
ylabel('average pressure [Pa]')
legend('matrix', 'fracture')
title('Results for the intrinsic fracture stiffness simulation')

figure
tplot = 1;
for i=1:nx
    X_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,1)';
    Y_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,2)';
    R_(i,:) = p_f(nx*(i-1)+1:nx*i,tplot);
end
pcolor(X_, Y_, R_);
title('fracture pressure')
shading interp
set(gca,'YDir','normal') 
colorbar

figure
for i=1:nx
    X_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,1)';
    Y_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,2)';
    R_2(i,:) = p_m(nx*(i-1)+1:nx*i,tplot);
end
pcolor(X_, Y_, R_2);
title('matrix pressure')
shading interp
set(gca,'YDir','normal') 
colorbar

figure
pcolor(X_, Y_, R_ - R_2);
title('Pf - Pm')
shading interp
set(gca,'YDir','normal') 
colorbar


figure
for i=1:nx
    X_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,1)';
    Y_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,2)';
    R_(i, :) = states{tplot}.stress(nx*(i-1)+1:nx*i,1);
end
pcolor(X_, Y_, R_);
shading interp
title('stress_x')
set(gca,'YDir','normal') 
colorbar

figure
for i=1:nx
    X_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,1)';
    Y_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,2)';
    R_(i, :) = states{tplot}.stress(nx*(i-1)+1:nx*i, 2);
end
pcolor(X_, Y_, R_);
shading interp
title('stress_y')
set(gca,'YDir','normal') 
colorbar


figure
for i=1:nx+1
    X__(i, :) =  G.nodes.coords((nx+1)*(i-1)+1:(nx+1)*i,1)';
    Y__(i, :) =  G.nodes.coords((nx+1)*(i-1)+1:(nx+1)*i,2)';
    R__(i,:) = u((nx+1)*2*(i-1)+1:2:(nx+1)*2*i,tplot);
end
pcolor(X__, Y__, R__);
shading interp
title('Ux')
set(gca,'YDir','normal') 
colorbar


figure
for i=1:nx+1
    X__(i, :) =  G.nodes.coords((nx+1)*(i-1)+1:(nx+1)*i,1)';
    Y__(i, :) =  G.nodes.coords((nx+1)*(i-1)+1:(nx+1)*i,2)';
    R__(i, :) = u((nx+1)*2*(i-1)+2:2:(nx+1)*2*i+1,tplot);
end
pcolor(X__, Y__, R__);
shading interp
title('Uy')
set(gca,'YDir','normal') 
colorbar


figure
pcolor(log10(perm_matrix));
title('log matrix perm')
shading interp
set(gca,'YDir','normal') 
colorbar

figure
pcolor(log10(perm_fracture));
title('log fracture perm')
shading interp
set(gca,'YDir','normal') 
colorbar
        
figure
pcolor(log10(reshape(E_f, nx,ny)));
title('Ef')
shading interp
set(gca,'YDir','normal') 
colorbar

figure
pcolor(log10(reshape(E_m,nx,ny)));
title('Em')
shading interp
set(gca,'YDir','normal') 
colorbar