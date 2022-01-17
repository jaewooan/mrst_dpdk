classdef TotalThermalEnergy < ComponentTotalMass
    methods
        function gp = TotalThermalEnergy(model, varargin)
            gp = gp@ComponentTotalMass(model, varargin{:});
            gp.dependencies = {};
            gp = gp.dependsOn({'PhaseThermalEnergy', 'RockThermalEnergy'});
            gp.label = 'E';
        end
        
        function energy = evaluateOnDomain(prop, model, state)
            [phaseEnergy, rockEnergy] ...
                = prop.getEvaluatedDependencies(state, 'PhaseThermalEnergy', ...
                                                       'RockThermalEnergy' );
            energy = rockEnergy;
            nph = model.getNumberOfPhases();
            for i = 1:nph
                energy = energy + phaseEnergy{i};
            end
            energy = prop.ensureMinimumDerivatives({energy});
            energy = energy{1};
        end
    end
end