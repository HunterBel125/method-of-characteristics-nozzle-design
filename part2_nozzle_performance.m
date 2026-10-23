% RECONSTRUCTED COMPONENT: Operating conditions and performance.
% Extracted from the completed integrated project.
%
% Requires part1_nozzle_geometry.m in the same folder.

clear; clc; close all

%% USER INPUTS

N_final = 50;

R = 287;

width_actual = 0.500;
exit_height_actual = 0.350;

psia_to_Pa = 6894.75729;

pb_psia = 0.1;
p0_max_psia = 500;

pb = pb_psia*psia_to_Pa;
p0_max = p0_max_psia*psia_to_Pa;

T_exit_min = 55;

%% GET GEOMETRY AND FLOW SOLUTION

result = part1_nozzle_geometry(N_final, false);

gamma = result.gamma;
M_exit = result.M_exit;
ht = result.ht;

x0 = result.x0;
y0 = result.y0;

x_wall = result.x_wall;
y_wall = result.y_wall;
M_wall = result.M_wall;

x_center = result.x_center;
M_center = result.M_center;

mesh_x = result.mesh_x;
mesh_y = result.mesh_y;
mesh_M = result.mesh_M;

L_triangle = x_wall(end) - x_center(end);

%% OPERATING CONDITIONS

Pcr3 = ...
    (1 + ((gamma-1)/2)*M_exit^2)^(-gamma/(gamma-1));

p0_min = pb/Pcr3;
p0_min_psia = p0_min/psia_to_Pa;

T0_min = T_exit_min*(1 + ((gamma-1)/2)*M_exit^2);
T0 = T0_min;

fprintf('================ OPERATING LIMITS ================\n')
fprintf('Back pressure, pb:                      %.6f psia\n', pb_psia)
fprintf('Minimum total pressure, p0_min:         %.6f psia\n', p0_min_psia)
fprintf('Maximum total pressure, p0_max:         %.6f psia\n', p0_max_psia)
fprintf('Minimum exit static temperature:        %.6f K\n', T_exit_min)
fprintf('Minimum required total temperature:     %.6f K\n\n', T0)

%% FULL WALL CONTOUR INCLUDING THROAT ARC

x_wall_full = [0, x0, x_wall(2:end)];
y_wall_full = [ht, y0, y_wall(2:end)];

[x_wall_unique, unique_idx] = unique(x_wall_full, 'stable');
y_wall_unique = y_wall_full(unique_idx);

%% INTERPOLATION GRID

nx = 350;
ny = 200;

xq = linspace(0, x_wall_full(end), nx);
yq = linspace(0, max(y_wall_full), ny);

[Xq, Yq] = meshgrid(xq, yq);

y_wall_interp = interp1( ...
    x_wall_unique, y_wall_unique, xq, 'linear', 'extrap');

Ywallq = repmat(y_wall_interp, ny, 1);

inside_nozzle = Yq >= 0 & Yq <= Ywallq;

%% CENTERLINE DISTRIBUTIONS

T_T0_center = ...
    1 ./ (1 + ((gamma-1)/2).*M_center.^2);

p_p0_center = T_T0_center.^(gamma/(gamma-1));

figure
plot(x_center, M_center, 'LineWidth', 2)
grid on
xlabel('x / h_t')
ylabel('Mach Number')
title('Centerline Mach Number Distribution')

figure
plot(x_center, p_p0_center, 'LineWidth', 2)
grid on
xlabel('x / h_t')
ylabel('p / p_0')
title('Centerline Pressure Ratio Distribution')

figure
plot(x_center, T_T0_center, 'LineWidth', 2)
grid on
xlabel('x / h_t')
ylabel('T / T_0')
title('Centerline Temperature Ratio Distribution')

%% TEMPERATURE CONTOUR MAP

T_mesh = T0 ./ ...
    (1 + ((gamma-1)/2).*mesh_M.^2);

F_T = scatteredInterpolant( ...
    mesh_x(:), mesh_y(:), T_mesh(:), 'natural', 'none');

Tq = F_T(Xq, Yq);
Tq(~inside_nozzle) = NaN;

figure
hold on

contourf(Xq, Yq, Tq, 40, 'LineColor', 'none')
colorbar
colormap turbo

plot(x_wall_full, y_wall_full, 'k-', 'LineWidth', 2)
plot([0 x_wall_full(end)], [0 0], 'k--', 'LineWidth', 1)

axis equal
grid on

xlabel('x / h_t')
ylabel('y / h_t')
title(sprintf( ...
    'Temperature Contour Map, T_0 = %.1f K, N = %d', ...
    T0, N_final))

%% PRESSURE CONTOUR MAP

p0_contour = p0_min;

p_mesh = p0_contour .* ...
    (1 + ((gamma-1)/2).*mesh_M.^2) ...
    .^(-gamma/(gamma-1));

F_p = scatteredInterpolant( ...
    mesh_x(:), mesh_y(:), p_mesh(:), 'natural', 'none');

pq = F_p(Xq, Yq);
pq(~inside_nozzle) = NaN;

figure
hold on

contourf(Xq, Yq, pq, 40, 'LineColor', 'none')
colorbar
colormap turbo

plot(x_wall_full, y_wall_full, 'k-', 'LineWidth', 2)
plot([0 x_wall_full(end)], [0 0], 'k--', 'LineWidth', 1)

axis equal
grid on

xlabel('x / h_t')
ylabel('y / h_t')
title(sprintf( ...
    'Pressure Contour Map, p_0 = %.1f psia, N = %d', ...
    p0_contour/psia_to_Pa, N_final))

%% PHYSICAL SCALING

scale_length = exit_height_actual/y_wall(end);

throat_height_actual = ht*scale_length;
nozzle_length_actual = x_wall(end)*scale_length;
triangle_length_actual = L_triangle*scale_length;

A_throat = width_actual*throat_height_actual;
A_exit = width_actual*exit_height_actual;

fprintf('================ DIMENSIONAL GEOMETRY ================\n')
fprintf('Specified exit height:                  %.6f m\n', exit_height_actual)
fprintf('Specified nozzle width:                 %.6f m\n', width_actual)
fprintf('Calculated throat height:               %.6f m\n', throat_height_actual)
fprintf('Calculated nozzle length:               %.6f m\n', nozzle_length_actual)
fprintf('Calculated uniform triangle length:     %.6f m\n\n', triangle_length_actual)

%% MASS FLOW RATE OVER OPERATING RANGE

p0_sweep_psia = linspace(p0_min_psia, p0_max_psia, 200);
p0_sweep = p0_sweep_psia*psia_to_Pa;

mdot_sweep = p0_sweep.*A_throat.*sqrt(gamma/(R*T0)) .* ...
    (2/(gamma+1))^((gamma+1)/(2*(gamma-1)));

figure
plot(p0_sweep_psia, mdot_sweep, 'LineWidth', 2)
grid on
xlabel('Total pressure, p_0 (psia)')
ylabel('Mass flow rate, mdot (kg/s)')
title('Mass Flow Rate Across Operating Pressure Range')

fprintf('================ MASS FLOW RATE RANGE ================\n')
fprintf('Mass flow at p0_min = %.3f psia: %.6f kg/s\n', ...
    p0_min_psia, mdot_sweep(1))
fprintf('Mass flow at p0_max = %.3f psia: %.6f kg/s\n\n', ...
    p0_max_psia, mdot_sweep(end))

%% FORCE ESTIMATES AT MAXIMUM PRESSURE

p0_force = p0_max;

T_exit = T0/(1 + ((gamma-1)/2)*M_exit^2);

p_exit = p0_force * ...
    (1 + ((gamma-1)/2)*M_exit^2)^(-gamma/(gamma-1));

a_exit = sqrt(gamma*R*T_exit);
V_exit = M_exit*a_exit;

mdot_force = p0_force*A_throat*sqrt(gamma/(R*T0)) * ...
    (2/(gamma+1))^((gamma+1)/(2*(gamma-1)));

T_throat = T0/(1 + ((gamma-1)/2)*1^2);
a_throat = sqrt(gamma*R*T_throat);
V_throat = a_throat;

% Original calculation retained.
thrust = mdot_force*(V_exit);

x_wall_dim = x_wall*scale_length;

M_wall_local = M_wall;

p_wall = p0_force .* ...
    (1 + ((gamma-1)/2).*M_wall_local.^2) ...
    .^(-gamma/(gamma-1));

F_vertical_signed = ...
    -width_actual*trapz(x_wall_dim, p_wall);

F_vertical_magnitude = abs(F_vertical_signed);

fprintf('================ FLOW / FORCE ESTIMATES AT 500 PSIA ================\n')
fprintf('Total temperature, T0:                  %.3f K\n', T0)
fprintf('Total pressure for force estimate:      %.3f psia\n', p0_force/psia_to_Pa)
fprintf('Back pressure:                          %.3f psia\n', pb/psia_to_Pa)
fprintf('Exit static temperature:                %.3f K\n', T_exit)
fprintf('Exit static pressure:                   %.3f psia\n', p_exit/psia_to_Pa)
fprintf('Exit velocity:                          %.3f m/s\n', V_exit)
fprintf('Throat velocity:                        %.3f m/s\n', V_throat)
fprintf('Mass flow rate at 500 psia:              %.6f kg/s\n', mdot_force)
fprintf('Nozzle thrust estimate:                 %.3f N\n', thrust)
fprintf('Vertical pressure force on one wall:     %.3f N downward\n\n', ...
    F_vertical_magnitude)