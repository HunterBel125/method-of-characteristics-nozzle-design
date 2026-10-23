% RECONSTRUCTED COMPONENT: Grid convergence study.
% Extracted from the completed integrated project.
%
% Requires part1_nozzle_geometry.m in the same folder.

clear; clc; close all

%% USER INPUTS

N_values = [5 10 15 20 30 40 50 60 80];

%% GET COMMON DESIGN PARAMETERS

% Obtain the design parameters from Part 1.
reference = part1_nozzle_geometry(N_values(1), false);

gamma = reference.gamma;
M_exit = reference.M_exit;
ht = reference.ht;

%% THEORETICAL AREA RATIO

Ae_At_expected = (1/M_exit) * ...
    ((2/(gamma+1)) * ...
    (1 + ((gamma-1)/2)*M_exit^2)) ...
    ^((gamma+1)/(2*(gamma-1)));

%% GRID CONVERGENCE STUDY

area_error = zeros(size(N_values));
exit_height = zeros(size(N_values));
nozzle_length = zeros(size(N_values));

for q = 1:length(N_values)

    N = N_values(q);

    result_q = part1_nozzle_geometry(N, false);

    area_error(q) = ...
        abs(result_q.y_wall(end)/ht - Ae_At_expected);

    exit_height(q) = result_q.y_wall(end);
    nozzle_length(q) = result_q.x_wall(end);
end

%% PRINT RESULTS

fprintf('================ GRID CONVERGENCE RESULTS ================\n')

fprintf('Target exit Mach number: %.6f\n', M_exit)
fprintf('Expected area ratio:     %.8f\n\n', Ae_At_expected)

fprintf('%8s %20s %20s %20s\n', ...
    'N', 'Area Ratio Error', 'Exit Height', 'Nozzle Length')

for q = 1:length(N_values)

    fprintf('%8d %20.8f %20.8f %20.8f\n', ...
        N_values(q), ...
        area_error(q), ...
        exit_height(q), ...
        nozzle_length(q))
end

fprintf('\n')

%% PLOT AREA-RATIO ERROR

figure

plot(N_values, area_error, 'o-', 'LineWidth', 2)

grid on
xlabel('N')
ylabel('|y_e/h_t - A_e/A_t|')
title('Grid Convergence: Area Ratio Error')