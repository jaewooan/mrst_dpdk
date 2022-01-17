classdef RockThermalEnergy < StateFunction

    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = RockThermalEnergy(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('RockInternalEnergy');
            gp = gp.dependsOn('RockMass', 'PVTPropertyFunctions');
            gp.label = 'E_R';
        end
        
        %-----------------------------------------------------------------%
        function uR = evaluateOnDomain(prop, model, state)
            uR    = prop.getEvaluatedDependencies(state, 'RockInternalEnergy');
            massR = model.getProp(state, 'RockMass');
            uR    = massR.*uR;
        end 
    end
    
end