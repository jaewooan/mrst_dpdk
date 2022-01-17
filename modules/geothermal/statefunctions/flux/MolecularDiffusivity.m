classdef MolecularDiffusivity < StateFunction
   
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function md = MolecularDiffusivity(model, varargin)
            md@StateFunction(model, varargin{:});
            md = md.dependsOn({'Temperature', 'PoreVolume'}, 'PVTPropertyFunctions');
            md = md.dependsOn({'pressure'}, 'state');
            md.label = 'd_i';
        end
        
        %-----------------------------------------------------------------%
        function d = evaluateOnDomain(prop, model, state)
            % Get dependencies
            [p, T, pv] = model.getProps(state, 'pressure', 'Temperature', 'PoreVolume');
            % Compute tourtuosity and porosity
            tau  = model.rock.tau;
            poro = pv./model.G.cells.volumes;
            % Comoute diffusivity
            ncomp = model.getNumberOfComponents();
            d     = cell(ncomp,1);
            for i = 1:ncomp
                if model.Components{i}.molecularDiffusivity > 0
                    d{i} = model.Components{i}.molecularDiffusivity.*tau.*poro;
                end
            end
        end
    end
    
end