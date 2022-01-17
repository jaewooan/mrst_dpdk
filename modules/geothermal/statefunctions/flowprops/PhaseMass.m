classdef PhaseMass < StateFunction
       
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = PhaseMass(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'Density', 'PoreVolume'}, 'PVTPropertyFunctions');
            gp = gp.dependsOn('s', 'state');
            gp.label = 'M_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function mass = evaluateOnDomain(prop, model, state)
            [rho, pv, s] = model.getProps(state, 'Density', 'PoreVolume', 's');
            nph  = model.getNumberOfPhases();
            mass = cell(1, nph);
            for i = 1:nph
                mass{i} = pv.*rho{i}.*s(:, i);
            end
        end       
    end
end
