classdef WellConductiveHeatFlux < StateFunction
   
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = WellConductiveHeatFlux(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp = gp.dependsOn({'Temperature'}, 'PVTPropertyFunctions');
            gp = gp.dependsOn({'RockHeatTransmissibility', 'FluidHeatTransmissibility'}, 'FlowDiscretization');
            gp.label = 'q_{h,c}';
        end
        
        %-----------------------------------------------------------------%
        function qh = evaluateOnDomain(prop, model, state)
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            [T, lambdaR, lambdaF] = model.ReservoirModel.getProps(state, 'Temperature', 'RockThermalConductivity', 'FluidThermalConductivity');
            T   = T(map.cells);
            isT = strcmpi(state.FacilityState.names, 'well_temperature');
            if any(isT)
                Tw = state.FacilityState.primaryVariables{isT};
            else
                Tw = vertcat(map.W.T);
            end
            Tw = Tw(map.perf2well);
            % Assume thermal conductivity tensors are diagonal for now
            G      = model.ReservoirModel.G;
            lambda = lambdaR + lambdaF;
            rock   = struct('perm', ones(numel(value(lambda)), 1));
            wi     = cell(numel(map.W),1);
            for i = 1:numel(map.W)
                radius = map.W(i).r;
                cells  = map.W(i).cells;
                wi{i}  = computeWellIndex(G, rock, radius, cells);
                wi{i}  = wi{i}.*lambda(cells);
            end
            wi = vertcat(wi{:});
            qh  = -wi.*(T - Tw);
        end
    end
    
end