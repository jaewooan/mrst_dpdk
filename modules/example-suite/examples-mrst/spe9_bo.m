function [description, options, state0, model, schedule, plotOptions] = spe9_bo(varargin)
%Example from the example suite, see description below.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`,

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    % One-line description
    description = ...
        ['SPE9 benchmark, see Killough, J.E., "Ninth SPE comparative  ' , ...
         'solution project: A reexamination of black-oil simulation.", ', ...
         'SPE Reservoir Simulation Symposium, 1995, SPE 29110-MS, '     , ...
         'doi: 10.2118/29110-MS'                                        ];
    % Optional input arguments
    options = struct();
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
    % Define module dependencies
    require ad-core ad-props ad-blackoil deckformat
    % Get initial state, model and schedule
    pth = getDatasetPath('spe9');
    fn  = fullfile(pth, 'BENCH_SPE9.DATA');
    [state0, model, schedule] = initEclipseProblemAD(fn, 'useMex', true, 'rowMajorAD', true);
    % Set plot options
    plotOptions = {'PlotBoxAspectRatio', [1,1,0.25], ...
                   'View'              , [120, 20] };
end