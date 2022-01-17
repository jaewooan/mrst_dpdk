classdef PhaseThermalEnergy < StateFunction
       
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = PhaseThermalEnergy(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'ComponentPhaseMass', 'PhaseInternalEnergy'});
            gp.label = 'E_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function energy = evaluateOnDomain(prop,model, state)
            [mass, u] = prop.getEvaluatedDependencies(state, 'ComponentPhaseMass' , ...
                                                             'PhaseInternalEnergy'); 
            nph    = model.getNumberOfPhases();
            ncomp  = model.getNumberOfComponents();
            energy = cell(1, nph);
            for i = 1:nph
                phmass = 0;
                for j = 1:ncomp
                    if isempty(mass{j,i}), continue; end
                    phmass = phmass + mass{j,i};
                end
                energy{i} = phmass.*u{i};
            end
        end

    end
end

