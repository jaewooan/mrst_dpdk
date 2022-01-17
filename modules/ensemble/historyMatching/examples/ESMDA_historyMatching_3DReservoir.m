%% History matching of 3D reservoir simulation
% This is a minimal example (and a sanity test) for creating an ensemble of
% 3D reservoirs with two phases (oil and water).
mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble 

mrstVerbose off

%% Set up and simulate the true solution
% We will here use an identical twin experiment, where we use the same
% problem for both generating the truth and as a base for our ensemble.

trueProblemName = 'ensemble_base_problem_3d_reservoir';
baseProblemOptions = {};


directoryTruth = fullfile(mrstOutputDirectory(), ...
                          'historyMatching', 'truth', ...
                          trueProblemName);
                      
topDirectory = fullfile(mrstOutputDirectory(), ...
                        'historyMatching', 'esmdaTutorial', trueProblemName);
                      
trueExample = MRSTExample(trueProblemName);
trueProblem = trueExample.getPackedSimulationProblem('Directory', directoryTruth);

plotExample = false;
rerunTrueProblemFromScratch = false;
overwriteObservation = true;

if rerunTrueProblemFromScratch
    clearPackedSimulatorOutput(trueProblem);
end
simulatePackedProblem(trueProblem);
if plotExample
    [wellSols, states, reports] = getPackedSimulatorOutput(trueProblem);
    trueExample.plot(states);
end

%% Generate observations
% Define a QoI object for storing the relevant observations we will use for
% history matching

trueQoI = WellQoIHM(...
    'wellNames', {'P1', 'P2'}, ...
    'names', {'qOs', 'qWs'}, ...
    'cumulative', false);

trueQoI = trueQoI.validateQoI(trueProblem);
trueObservations = trueQoI.getQoI(trueProblem);

% Define observation uncertainty 
obsStdDev = 0.0004; %*0.1;

% Create a separate ResultHandler for the observations 
observationResultHandler = trueQoI.ResultHandler;
observationResultHandler.dataPrefix = 'observedQoI';

% Add some observation noise and store output
if numel(observationResultHandler.getValidIds) < 1 || overwriteObservation
    for w = 1:numel(trueQoI.wellNames)
        perturbedObservations(w) = trueObservations(w);
        for f = 1:numel(trueQoI.names)
            trueVals = trueObservations(w).(trueQoI.names{f});
            perturbedObservations(w).(trueQoI.names{f}) = trueVals + randn(size(trueVals))*obsStdDev;
        end
    end
    observationResultHandler{1} = {perturbedObservations};
end



%% Select and populate samples for stochastic configurations class

ensembleSize = 70;


rockData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    rockData{i}.poro = gaussianField(trueExample.model.G.cartDims, [0.2 0.4]); 
    rockData{i}.perm = rockData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-rockData{i}.poro).^2);
end

rockSamples = RockSamplesHM('data', rockData);

%% Select quantity of interest class matching the what we have as observations
% We validate the QoI with the trueProblem, since this will be our ensemble
% base problem as well.

qoi = WellQoIHM('wellNames', {'P1', 'P2'}, ...
              'names', {'qOs', 'qWs'}, ...
              'observationResultHandler', observationResultHandler, ...
              'observationCov', obsStdDev^2);


%% Create the ensemble
esmdaEnsemble = MRSTHistoryMatchingEnsemble(trueExample, rockSamples, qoi, ...
    'alpha', [28/3 7 4 2], ...
    'directory', fullfile(topDirectory, 'esmda'), ...
    'simulationStrategy', 'spmd', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true)

%% Displaying the observations and observation error cov through the ensemble
disp('observation and scaling vector')
[obsVector, obsScaling] = esmdaEnsemble.qoi.getObservationAndScaling()
disp('observation error covariance matrix')
esmdaEnsemble.qoi.getObservationErrorCov()

%% Run ensemble
esmdaEnsemble.simulateEnsembleMembers();

%% Get simulated observations
disp('simulated observations')
size(esmdaEnsemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
size(esmdaEnsemble.getEnsembleSamples())

%% Do history matching
disp('updated sample object:')
esmdaEnsemble.doHistoryMatching();

%% Run resulting forecast
esmdaEnsemble.simulateEnsembleMembers();


%% Plot results
esmdaEnsemble.plotQoI('subplots', true, ...
                      'clearFigure', false, ...
                      'subIterations', true, ...
                      'legend', {'observations', 'posterior', 'ES-MDA it 3',...
                                 'ES-MDA it 2', 'ES-MDA it 1', 'prior'});







%% EnKF
%%%%%%%%%%%%%%%%%%%%%

enkfEnsemble = MRSTHistoryMatchingEnsemble(trueExample, rockSamples, qoi, ...
    'alpha', 1, ... %'alpha', [2 2] ,... %[28/3 7 4 2], ...
    'directory', fullfile(topDirectory, 'enkf'), ...
    'simulationStrategy', 'spmd', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true, ...
    'verboseSimulation', false)


%enkfEnsemble.simulateEnsembleMembers('range', 1);


%% Loop 

totalNumTimesteps = 15; %numel(trueProblem.SimulatorSetup.schedule.step.val);
enkfInvervals = {(1:2), (3:4), (5:6), (7:8), (9:10)};

for i = 1:numel(enkfInvervals)
    enkfEnsemble.updateHistoryMatchingInterval(enkfInvervals{i});
    enkfEnsemble.simulateEnsembleMembers();
    enkfEnsemble.doHistoryMatching();
end


%% Forecast the posterior
enkfEnsemble.updateHistoryMatchingInterval(1:15);
enkfEnsemble.simulateEnsembleMembers();



%% Plot history matching
enkfEnsemble.plotQoI('subplots', true, ...
                     'clearFigure', false, ...
                     'subIterations', true, ...
                     'observationIndices', (1:10), ...
                     'legend', {'observations', 'posterior', 'EnKF it 4',...
                                 'EnKF it 3', 'EnKF it 2', 'EnKF it 1', ...
                                 'prior'});

%% Copyright Notice
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
