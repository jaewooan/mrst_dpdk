clear
clc
close all
addpath('./SGS/')

%% Input from simulation: please check if this dimension is consistent with the values in the simulation!
Lx = 100; Ly = 100;
nx = 10; ny = 10; % check 
nCase = 2; % number of realization

%% SGS Setting
covar.model = 'gaussian'; % type of covariance function, see functions/kriginginitiaite.m for option
covar.range0 = [nx/10 nx/10]; % range of covariance [d_nx d_ny]
covar.azimuth = 0; % orientation of the covariance
covar.c0 = 1; % variance
covar.alpha = 1; % parameter of covariance function (facult)
parm.saveit = false; % save in a file
parm.seed_path = 'shuffle'; % seed for the path
parm.seed_search = 'shuffle'; % seed for the search (equal distance node)
parm.seed_U = 'shuffle'; % seed of the node value sampling
parm.mg = 1 ; % multigrid vs randomized path
neigh.wradius = 3; % maximum range of neighboorhood search, factor of covar.range.
neigh.lookup = false; % lookup table.
neigh.nb = 40; % maximum number of neighborhood
method = 'trad'; % SGS method: cst_par, cst_par_cond, hybrid, varcovar

%% SGS
[logPermM, t1] = SGS(nx, ny, nCase, covar, neigh, parm, method);
[logPermF, t2] = SGS(nx, ny, nCase, covar, neigh, parm, method);
[logEM, t1] = SGS(nx, ny, nCase, covar, neigh, parm, method);
[logEF, t2] = SGS(nx, ny, nCase, covar, neigh, parm, method);

i_index = rem([1:nx*ny]'-1, nx) + 1;
j_index = floor(([1:nx*ny]'-1)/nx)+1;
x_centroid = Lx/nx*(i_index-1/2);
y_centroid = Ly/ny*(j_index-1/2);
perm_m_tot = exp(reshape(logPermM, nx*ny, nCase));
perm_f_tot = 100*exp(reshape(logPermF, nx*ny, nCase));
EM_tot = 1e9*exp(reshape(logEM, nx*ny, nCase));
EF_tot = 1e7*exp(reshape(logEF, nx*ny, nCase));

%%
save('perm.mat', 'perm_m_tot', 'perm_f_tot', 'EM_tot', 'EF_tot', 'i_index', 'j_index', 'x_centroid', 'y_centroid', 'nCase');

figure
perm1 = exp(logPermM(:,:,1));
pcolor(logPermM(:,:,1))
shading interp
colorbar

figure
perm2 = 100*exp(logPermF(:,:,1));
pcolor(logPermF(:,:,1))
shading interp
colorbar