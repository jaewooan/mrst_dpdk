function [p_m, p_f, u, states] = simDC_mech(state0, dt, model, bc_f0)
%
% SYNOPSIS:
%   function [p_m, p_f, u, states] = simDC_mech(states, dt, end_time, model, bc_f0)
%
% DESCRIPTION:
%    A function which runs a consolidation type simulation for the
%    dual-continuum (DC) poroelastic model. Includes the induced problem  
%    at t = 0, in which the applied stress induces a pressure within the
%    DC domain. 
%
% PARAMETERS:
%   states     - struc of states of the system
%   dt         - time step
%   end_time   - end simulation time
%   model      - instance of a dual-perm mechanical class
%   bc_f0      - fluid boundary conditions for induced problem
%
% RETURNS:
%   p_m (matrix pressure), p_f (fracture pressure), u (displacements),
%   states
%
% SEE ALSO: example_void_fractures
%
    
    %% Solver and simulation parameters
    solver = NonLinearSolver();
    n = 1 + length(dt); % additional 1 because we consider the induced problem
                         % to be at t = 0. 
    p_m = zeros(model.G.cells.num, n);
    p_f = zeros(model.G.cells.num, n);
    u = zeros(model.G.nodes.num*model.G.griddim, n);
    states = cell(n, 1);   
    
    W = [];
    W = addWell(W, model.G, model.rock, [5,15,25,35, 45, 55,65,75], ...
        'comp_i', [0.17, 0.74, 0.09],'Name', ['Prod1'], 'Val',...
        1, 'sign', -1, 'Type', 'bhp','Radius', 0.01,'Dir','x');

    state0.wellSol = initWellSolAD([], model, state0);
    state0.xd = zeros(nnz(~model.mechModel.operators.isdirdofs), 1);
    state0 = addDerivedQuantities(model.mechModel, state0);

    % define facilities model, which for this case is empty
    model.FacilityModel = FacilityModel(model.fluidModel); % map facility model from fluidModel up to main model
    model.FacilityModel = model.FacilityModel.setupWells(W);
    
    %% Loop and solve
    for i = 1:n
        if i == 1
        % induced problem
            state = solver.solveTimestep(state0, dt(1), model, 'bc', bc_f0);
            states{i} = state;
            p_m(:,i) = state.pressure_matrix;
            p_f(:,i) = state.pressure;
            u(:,i) = state.u;
            bc_f = fluxside([], model.G, 'WEST', 0, 'sat', 1);
            bc_f = fluxside(bc_f, model.G, 'EAST', 0, 'sat', 1);
            bc_f = fluxside(bc_f, model.G, 'SOUTH', 0, 'sat', 1);
            bc_f = pside(bc_f, model.G, 'NORTH', 0, 'sat', 1);
        else
            state = solver.solveTimestep(states{i-1}, dt(i-1), model, 'bc', bc_f);
            states{i} = state;
            p_m(:,i) = state.pressure_matrix;
            p_f(:,i) = state.pressure;
            u(:,i) = state.u;
        end
    end
    
end