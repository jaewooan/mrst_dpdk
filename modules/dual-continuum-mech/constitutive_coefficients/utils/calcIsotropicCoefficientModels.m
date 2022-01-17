function [B_f, B_m, invN_f, invQ, invN_m] = calcIsotropicCoefficientModels(coefficients, model)
%
% SYNOPSIS:
%   function [B_f, B_m, invN_f, invQ, invN_m] = isotropicCoefficientModels(coefficients, model)
%
% DESCRIPTION: 
%   Computes anisotropic constitutive coefficients necessary for the
%   poroelastic dual-constitutive constitutive model using
%   the coefficient models presented in Ashworth and Doster 
%   (TBR)
%
% PARAMETERS:
%   coefficients - class instance
%   model 
%
% RETURNS:
%   B_f    - Effective Biot coefficient of the fracture
%   B_f    - Effective Biot coefficient of the matrix
%   invN_f - Effective solid Biot modulus of the fracture
%   invQ   - Inter-continuum pressure cross coupling coefficient
%   invN_m - Effective solid Biot modulus of the matrix
%
% EXAMPLE:
%
% SEE ALSO: isotropicStiffCoefficients, ConstitutiveCoefficients
%           
%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

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
C = model.mech.C;
d = model.G.griddim;
if d == 3
    rank2 = @(x) bsxfun(@times, x, [1, 1, 1, 0, 0, 0]);
else
    rank2 = @(x) bsxfun(@times, x, [1, 1, 0]);
end
identity = rank2(ones(model.G.cells.num, 1));

K_dr = (1/(d^2))*doubledot(identity, doubledot(identity, C, model.G), model.G); 
[v_m, v_f] = deal(model.rock_matrix.vol_fraction, model.rock.vol_fraction);
[K_m, b_m, invn_m, K_f, b_f, invn_f] = deal(coefficients.iP.C_m, coefficients.iP.b_m, ... 
                                            coefficients.iP.invn_m, coefficients.iP.C_f, ...
                                            coefficients.iP.b_f, coefficients.iP.invn_f);
% original equation
B_f_ = b_f - b_f.*((K_dr-K_f)./(K_m-K_f)); 
B_m_ = b_m.*((K_dr-K_f)./(K_m-K_f));
invN_f_ = (b_f./(K_m-K_f)).*(B_f - v_f.*b_f) + v_f.*invn_f;
invN_m_ = (b_m./(K_m-K_f)).*(v_m.*b_m - B_m) + v_m.*invn_m;
invQ_ = (b_m./(K_m-K_f)).*(v_f.*b_f - B_f);

% Jihoon Kim's formuilation (2012)
B_f = -b_f.*K_dr.*v_f./K_f;
B_m = -b_m.*K_dr.*v_m./K_m;
invn_Mm; %inverse of matrix biot modulus: 1/Mm = phi*cf +(b_m-phi)/Ks = phi*cf + invn_m
invn_Mf; %inverse of fracture biot modulus: 1/Mf = phi*cf +(b_f-phi)/Ks = phi*cf + invn_f
dff = v_f./K_f.*b_f./B
invN_f = (b_f./(K_m-K_f)).*(B_f - v_f.*b_f) + v_f.*invn_f;
invN_m = (b_m./(K_m-K_f)).*(v_m.*b_m - B_m) + v_m.*invn_m;
invQ = (b_m./(K_m-K_f)).*(v_f.*b_f - B_f);


B_m = rank2(B_m);
B_f = rank2(B_f);
end

