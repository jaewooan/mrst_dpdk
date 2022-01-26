classdef MandelAnalyticalSolution
    %MESHSTRUCTURED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type = 'mandel_analytical_solution_one_quarter';
        a
        b
        young_modulus
        poisson_ratio
        biot_coefficient
        biot_modulus
        permeability
        fluid_viscosity
        %
        force
        %
        drained_bulk_modulus
        shear_modulus
        uniaxial_vertical_compressibility
        undrained_bulk_modulus
        Skempton_coefficient
        undrained_poisson_ratio
        fluid_diffusivity
        characteristic_time
        %
        initial_pressure
        max_initial_displacement_X
        min_initial_displacement_Z
        max_final_displacement_X
        min_final_displacement_Z
        A1
        A2
        %
        time
        displacement
        pressure
    end
    
    methods
        function MaAnSo = MandelAnalyticalSolution(a, b, young_modulus, poisson_ratio, biot_coefficient, biot_modulus, permeability, fluid_viscosity,force)
            if nargin >= 1
                MaAnSo.a = a;
                MaAnSo.b = b;
                MaAnSo.young_modulus = young_modulus;
                MaAnSo.poisson_ratio = poisson_ratio;
                MaAnSo.biot_coefficient = biot_coefficient;
                MaAnSo.biot_modulus = biot_modulus;
                MaAnSo.permeability = permeability;
                MaAnSo.fluid_viscosity = fluid_viscosity;
                MaAnSo.force = force*a;
                %
                MaAnSo.drained_bulk_modulus = young_modulus / 3 / (1 - 2 * poisson_ratio);
                MaAnSo.shear_modulus = 0.5 * young_modulus / (1 + poisson_ratio);
                MaAnSo.uniaxial_vertical_compressibility = (1 + poisson_ratio) * (1 - 2 * poisson_ratio) / young_modulus / (1 - poisson_ratio);
                MaAnSo.undrained_bulk_modulus = MaAnSo.drained_bulk_modulus + biot_coefficient^2*biot_modulus;
                MaAnSo.Skempton_coefficient = biot_coefficient*biot_modulus / MaAnSo.undrained_bulk_modulus;
                MaAnSo.undrained_poisson_ratio = ( 3 * poisson_ratio + biot_coefficient * MaAnSo.Skempton_coefficient * (1 - 2 * poisson_ratio) ) / ( 3 - biot_coefficient * MaAnSo.Skempton_coefficient * (1 - 2 * poisson_ratio));
                MaAnSo.fluid_diffusivity = permeability / fluid_viscosity / (1 / biot_modulus + biot_coefficient^2 * MaAnSo.uniaxial_vertical_compressibility);
                MaAnSo.characteristic_time = a^2 / MaAnSo.fluid_diffusivity;
                %
                MaAnSo.initial_pressure = MaAnSo.Skempton_coefficient*(1 + MaAnSo.undrained_poisson_ratio) * abs(MaAnSo.force) / 3 / a;
                MaAnSo.max_initial_displacement_X = abs(MaAnSo.force) * MaAnSo.undrained_poisson_ratio / 2 / MaAnSo.shear_modulus;
                MaAnSo.min_initial_displacement_Z = -abs(MaAnSo.force) * (1 - MaAnSo.undrained_poisson_ratio) * b / 2 / MaAnSo.shear_modulus / a;
                MaAnSo.max_final_displacement_X = abs(MaAnSo.force) * poisson_ratio / 2 / MaAnSo.shear_modulus;
                MaAnSo.min_final_displacement_Z = -abs(MaAnSo.force) * (1 - poisson_ratio) * b / 2 / MaAnSo.shear_modulus / a;
                MaAnSo.A1 = 3 / (MaAnSo.Skempton_coefficient * ( 1 + MaAnSo.undrained_poisson_ratio) );
                MaAnSo.A2 = MaAnSo.biot_coefficient * (1 - 2 * MaAnSo.poisson_ratio) / (1 - MaAnSo.poisson_ratio);
            end
        end % MandelAnalyticalSolution         

        function MaAnSo = compute_pressure_solution(MaAnSo, time, pressure_coordinates)
            
            number_pressure_coordinates = size(pressure_coordinates,2);
            
            MaAnSo.pressure = zeros(number_pressure_coordinates, 1);
            
            if ( abs(time) <= eps)
                MaAnSo.pressure = MaAnSo.initial_pressure * ones(max(size(pressure_coordinates,2)), 1);
            else

                % Series to compute  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % sum1 = \sum_{n=1}^\infty 
                %           \frac{sin \alpha_n}{\alpha_n - sin alpha_n*cos\alpha_n}*

                dsum1 = @(an, x, a, cf, t) ...
                    sin(an) * (cos(an*x/a) - cos(an)) * exp(-an^2*cf*t/a^2) / ...
                    (an - sin(an)*cos(an));
               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                tol = eps;
                itmax = 400;

                for iCoord = 1:number_pressure_coordinates

                    %
                    % Compute analitical pressure at coordinate x and time t
                    x = pressure_coordinates(1,iCoord);

                    ds1 = 1000*eps;
                    sum1 = 0.0;
                    iter = 0;

                    while ( (abs(ds1) > tol) && (iter < itmax) )
                        iter = iter + 1;
                        %
                        % Compute coefficient alfa_n
                        resan = 100*tol;
                        iter_an = 0;
                        an_k = -pi/2 + iter*pi-100*eps;

                        while ( resan > tol && iter_an < itmax)
                            iter_an = iter_an + 1;
                            fx = (tan(an_k) - MaAnSo.A1/MaAnSo.A2*an_k);
                            dfx = (1 + tan(an_k)^2 - MaAnSo.A1/MaAnSo.A2);
                            an = an_k -fx/dfx;
                            resan = abs(an - an_k);
                            an_k = an;
                        end

                        ds1 = dsum1(an, x, MaAnSo.a, MaAnSo.fluid_diffusivity, time);
                        sum1 = sum1 + ds1;

                    end

                    MaAnSo.pressure(iCoord) = 2 * MaAnSo.initial_pressure * sum1;

                end
                
            end
        end % compute_solution
		
        function MaAnSo = compute_displacement_solution(MaAnSo, time, displacement_coordinates)
            
            number_displacement_coordinates = size(displacement_coordinates,2);
            
            MaAnSo.displacement = zeros(2 * number_displacement_coordinates, 1);
            
            if ( abs(time) <= eps)
                MaAnSo.displacement(1:2:end) = MaAnSo.max_initial_displacement_X / MaAnSo.a * displacement_coordinates(1,:)';
                MaAnSo.displacement(2:2:end) = MaAnSo.min_initial_displacement_Z / MaAnSo.b * displacement_coordinates(2,:)';
            else

                % Series to compute  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % sum1 = \sum_{n=1}^\infty 
                %           \frac{sin \alpha_n}{\alpha_n - sin alpha_n*cos\alpha_n}*

                dsum2 = @(an, a, cf, t) ...
                    sin(an) * cos(an) * exp(-an^2*cf*t/a^2) / (an - sin(an)*cos(an));

                dsum3 = @(an, x, a, cf, t) ...
                    cos(an) * sin(an*x/a) * exp(-an^2*cf*t/a^2) / (an - sin(an)*cos(an));                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                tol = eps;
                itmax = 400;

                for iCoord = 1:number_displacement_coordinates

                    %
                    % Compute analitical pressure at coordinate x and time t
                    x = displacement_coordinates(1,iCoord);
                    y = displacement_coordinates(2,iCoord);

                    ds2 = 100*eps;
                    ds3 = 100*eps;
                    sum2 = 0.0;
                    sum3 = 0.0;
                    iter = 0;

                    while ( ( (abs(ds2) > tol) || (abs(ds3) > tol) ) && (iter < itmax) )

                        iter = iter + 1;

                        %
                        % Compute coefficient alfa_n
                        resan = 100*tol;
                        iter_an = 0;
                        an_k = -pi/2 + iter*pi-100*eps;

                        while ( resan > tol && iter_an < itmax)
                            iter_an = iter_an + 1;
                            fx = (tan(an_k) - MaAnSo.A1/MaAnSo.A2*an_k);
                            dfx = (1 + tan(an_k)^2 - MaAnSo.A1/MaAnSo.A2);
                            an = an_k -fx/dfx;
                            resan = abs(an - an_k);
                            an_k = an;
                        end

                        if (abs(ds2) > tol)
                            ds2 = dsum2(an, MaAnSo.a, MaAnSo.fluid_diffusivity, time);
                            sum2 = sum2 + ds2;
                        end
                        if (abs(ds3) > tol)
                            ds3 = dsum3(an, x, MaAnSo.a, MaAnSo.fluid_diffusivity, time);
                            sum3 = sum3 + ds3;
                        end        

                    end
                    MaAnSo.displacement(iCoord*2 - 1) = ( abs(MaAnSo.force) / MaAnSo.shear_modulus) * ( (0.5 * MaAnSo.poisson_ratio - MaAnSo.undrained_poisson_ratio * sum2) * x / MaAnSo.a + sum3 );
                    MaAnSo.displacement(iCoord*2) = ( abs(MaAnSo.force) / MaAnSo.shear_modulus / MaAnSo.a) * ( -0.5 * (1 - MaAnSo.poisson_ratio) + (1 - MaAnSo.undrained_poisson_ratio) * sum2 ) * y;
                end
                
            end
        end % compute_solution        
    end
    
end
