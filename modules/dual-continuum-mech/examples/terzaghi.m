% % Copyright Notice
% 
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>



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
nx = 1; ny = 10; 
Lx = 1; Ly = 10; %m
km = 100; kf = 100; %md
porom = 0.5; porof = 0.5;
vol_f = 0.5; % frac vol fraction
mu = 1/centi/poise;%0.0981; % centipoise = 1e-3 pa
rho = 1;%1000; %kg/m3
cf = 1e-10;%4.4e-10; %fluid compressibility 1/pa
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

alpha = 1;
h = Ly/ny;         %discretized, m
K1D = K_m+4/3*G_shear;    %uniaxial drained bulk modulus
c_ = kf/nu_m*M_m * K1D/(K1D+alpha^2*M_m); %consolidation coefficient
tc = Ly^2/c_;       %characteristic time
Nt = 100;         %number of timesteps
N = 500;          %number of elements in series
m = [0:1:N];      %fourier series sequence
p0 = alpha*M_m*force/(K1D+alpha^2*M_m);  % undrained pressure, Pa
u0 = Ly*force/(K1D+alpha^2*M_m);        %max initial dispacement, m

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
              'poro', ones(G.cells.num, 1)*porof,... 
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
    if i == 3
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
bc_el_sides{2}.el_bc.disp_bc.mask(:, 2) = false;   % east x free, y free
bc_el_sides{3} = bc{3}; 
bc_el_sides{3}.el_bc.disp_bc.mask(:, :) = true;    % north x fixed, y fixed
bc_el_sides{4} = bc{4}; 
bc_el_sides{4}.el_bc.disp_bc.mask(:, 2) = false;   % south x free, y fixed

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
d1 = sqrt(2)*pi;%0.1; % spacing of fracture set 1
d2 = d1; % spacing of fracture set 2
fracture_spacing = repmat([d1,d2],G.cells.num,1);
shape_factor_name = 'Lim_AzizShapeFactor';
DC_model.transfer_model_object = SimpleTransferFunction(shape_factor_name, fracture_spacing);
DC_model = DC_model.validateModel();

%% Setup initial state and fluid BCs
pressure = ones(G.cells.num,1)*p0; % pa
u_init = zeros(2*(nx+1)*(ny+1), 1);
u_init(1:2:end) = 0;
u_init(2:2:end) = u0*(Ly-G.nodes.coords(:,2));

state0 = struct('pressure', pressure, 'pressure_matrix', pressure, 's', ones(G.cells.num, 1), 'swm', ones(G.cells.num, 1));
% state0 = struct('pressure', ones(G.cells.num,1)*p0, 'pressure_matrix', ones(G.cells.num,1)*p0, 's', ones(G.cells.num, 1), 'swm', ones(G.cells.num, 1),...
%                 'u', u_init);

% need to initiate the fluid bc's, bc's are the same for micro and macro scales 
%bc_f0 = pside([], G, 'WEST', 10000000, 'sat', 1);
bc_f0 = fluxside([], G, 'WEST', 0, 'sat', 1);
bc_f0 = fluxside(bc_f0, G, 'EAST', 0,'sat', 1);
bc_f0 = fluxside(bc_f0 , G, 'SOUTH', 0, 'sat', 1);
bc_f0 = pside(bc_f0, G, 'NORTH', 0, 'sat', 1);
%bc_f0 = fluxside(bc_f0, G, 'NORTH', 0, 'sat', 1);
%bc_f0 = pside(bc_f0 , G, 'SOUTH', 0, 'sat', 1);
%bc_f0 = pside(bc_f0 , G, 'NORTH', 0, 'sat', 1);


%% Simulate 
time = [1:1:Nt]*tc/Nt*100;
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
% close all
% close all
% figure
% for i=1:ny
%     X_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,1)';
%     Y_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,2)';
%     R_(i,:) = p_f(nx*(i-1)+1:nx*i,size(time,2));
% end
% pcolor(X_, Y_, R_);
% shading interp
% set(gca,'YDir','normal') 
% colorbar
% 
% figure
% for i=1:ny
%     X_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,1)';
%     Y_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,2)';
%     R_(i,:) = p_m(nx*(i-1)+1:nx*i,size(time,2));
% end
% pcolor(X_, Y_, R_);
% shading interp
% set(gca,'YDir','normal') 
% colorbar
% 
% 
% figure
% for i=1:ny
%     X_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,1)';
%     Y_(i, :) =  G.cells.centroids(nx*(i-1)+1:nx*i,2)';
%     R_(i, :) = states{1}.stress(nx*(i-1)+1:nx*i, 2);
% end
% pcolor(X_, Y_, R_);
% shading interp
% set(gca,'YDir','normal') 
% colorbar
% %%
% figure
% for i=1:ny+1
%     X__(i, :) =  G.nodes.coords((nx+1)*(i-1)+1:(nx+1)*i,1)';
%     Y__(i, :) =  G.nodes.coords((nx+1)*(i-1)+1:(nx+1)*i,2)';
%     R__(i,:) = u((nx+1)*2*(i-1)+1:2:(nx+1)*2*i,size(time,2));
% end
% pcolor(X__, Y__, R__);
% set(gca,'YDir','normal') 
% colorbar
% 
% 
% figure
% for i=1:ny+1
%     X__(i, :) =  G.nodes.coords((nx+1)*(i-1)+1:(nx+1)*i,1)';
%     Y__(i, :) =  G.nodes.coords((nx+1)*(i-1)+1:(nx+1)*i,2)';
%     R__(i, :) = u((nx+1)*2*(i-1)+2:2:(nx+1)*2*i+1,5);
% end
% pcolor(X__, Y__, R__);
% set(gca,'YDir','normal') 
% colorbar



%%
% Terzaghi: simple analytical solution (1D problem)
% author: Ruslan Rin (iskhakov@stanford.edu)
% 5/13/2016
clc
z = [h/2:h:Ly];    % center of our blocks
z2 = [0:h:Ly];
Nb = 5;
if (Nb < 20)
  %we need to shift since numerical schemes uses 9 blocks to mimic
  %boundaries capabilities
  z(Nb) = z(Nb-1);    
  z2(1) = -h/8; 
end;

%analytical solution
for i = 1:Nt  
  tmp_exp = exp(-pi^2*c_*time(i)*(2*m+1).^2/(4*Ly^2));
  for iz = 1:Nb
    xtmp = sum(1./(2*m+1).*tmp_exp.*sin((2*m+1)*pi.*z(iz)/2/Ly));
    p(iz,i) = 4/pi*p0*xtmp; %Pa
    xtmp = sum(1./(2*m+1).^2.*tmp_exp.*cos((2*m+1)*pi.*z2(iz)/2/Ly));
    u_a(iz,i) = alpha*p0/K1D*(Ly-z2(iz)-8*Ly/pi^2*xtmp) + u0;
  end;
end;
tt = time/tc;

% pressure and dispacement 
pmax = 85;
umax = 0.8;

figure; set(gcf,'color','white'); set(gcf, 'Position', [400 200 600 500])
plot(tt, u_a(1,:));
hold on
plot(tt, u(1,:));
ylabel('dispacement, m');
xlabel('time, t/t_c');
legend('u_{top,analytical}','u_{top,numerical}',...
          'Location', 'SouthWest');

figure;
plot(tt,p(Nb,:));
hold on
plot(tt,p_m(Nb,:));
ylabel('pressure, pa');
xlabel('time, t/t_c');
legend('p_{bot,analytical}', 'p_{bot,numerical}',...
          'Location', 'SouthWest');
%Image = getframe(gcf); imwrite(Image.cdata, 'terzagi10_comparison.png');




        