function gdz = setup_gdz(model, grad_op)
%Undocumented Utility Function

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

    g = model.getGravityVector();
    griddim = model.G.griddim;
    grad = zeros(sum(model.operators.internalConn), griddim);

    % FIX move loop
    for k = 1:griddim
        grad(:, k) = grad_op(model.G.cells.centroids(:, k));
    end

    gdz = grad * g';
end
