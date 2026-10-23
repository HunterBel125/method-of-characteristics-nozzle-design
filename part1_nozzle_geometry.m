function result = part1_nozzle_geometry(N, make_plots)
% PART1_NOZZLE_GEOMETRY
% Design a symmetric two-dimensional nozzle using the method of
% characteristics.
%
% The solver computes the upper half of the nozzle. The contour plots
% reflect this solution across the centerline to display the full nozzle.
%
% Examples:
%   result = part1_nozzle_geometry;
%   result = part1_nozzle_geometry(50, true);
%   result = part1_nozzle_geometry(50, false);
%
% Inputs:
%   N          - Number of initial characteristic points.
%   make_plots - true: display figures.
%                false: return results without displaying figures.
%
% Figures:
%   1. Upper-half geometry and characteristic network.
%   2. Full symmetric Mach contour with axes and colorbar.
%   3. Clean full symmetric Mach contour with outline and colorbar.
%
% Nothing is saved automatically.

%% Inputs and design conditions

if nargin < 1 || isempty(N)
    N = 50;
end

if nargin < 2 || isempty(make_plots)
    make_plots = true;
end

validateattributes(N, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2});

gamma  = 1.4;
M_exit = 6;

% Geometry is normalized by the throat half-height.
ht = 1;
Rt = ht;

% Contour resolution for the computed upper half.
n_colors = 60;
plot_nx  = 700;
plot_ny  = 250;

% Appearance.
wall_color = [0.07, 0.12, 0.20];
wall_line_width = 2.2;
end_line_width  = 1.3;

%% Solve the upper-half characteristic network

result = run_moc_nozzle(N, gamma, M_exit, ht, Rt);

result.N      = N;
result.gamma  = gamma;
result.M_exit = M_exit;
result.ht     = ht;
result.Rt     = Rt;

% Allow other project files to call the solver without plotting.
if ~make_plots
    return;
end

%% Extract results

x0 = result.x0;
y0 = result.y0;

x_wall     = result.x_wall;
y_wall     = result.y_wall;
theta_wall = result.theta_wall;
M_wall     = result.M_wall;

x_center = result.x_center;
M_center = result.M_center;

mesh_x = result.mesh_x;
mesh_y = result.mesh_y;
mesh_M = result.mesh_M;

segments = result.segments;

fprintf('\nMETHOD OF CHARACTERISTICS NOZZLE\n');
fprintf('--------------------------------\n');
fprintf('Design exit Mach number: %.2f\n', M_exit);
fprintf('Number of initial characteristics: %d\n', N);
fprintf('Exit half-height / throat half-height: %.6f\n', ...
    y_wall(end) / ht);
fprintf('Full exit height / throat half-height: %.6f\n', ...
    2 * y_wall(end) / ht);
fprintf('Nozzle length / throat half-height: %.6f\n', ...
    x_wall(end) / ht);

%% Construct the complete upper wall

% Include the throat, initial curved wall, and downstream designed wall.
x_wall_full = [0, x0, x_wall(2:end)];
y_wall_full = [ht, y0, y_wall(2:end)];

% Remove repeated x coordinates before interpolating the wall.
[x_wall_full, unique_indices] = unique(x_wall_full, 'stable');
y_wall_full = y_wall_full(unique_indices);

%% Identify the uniform exit-flow region

% The final characteristic connects the last centerline point to the
% exit wall. Under the assumed uniform axial exit-flow condition, the
% downstream triangle is assigned the design Mach number.
%
% Triangle vertices:
%   1. Final centerline point.
%   2. Exit wall point.
%   3. Exit-plane centerline point.

uniform_x = [ ...
    x_center(end), ...
    x_wall(end), ...
    x_wall(end)];

uniform_y = [ ...
    0, ...
    y_wall(end), ...
    0];

uniform_flow_supported = ...
    abs(M_center(end) - M_exit) < 1e-6 && ...
    abs(M_wall(end) - M_exit) < 1e-6 && ...
    abs(theta_wall(end)) < 1e-6 && ...
    x_wall(end) > x_center(end);

if ~uniform_flow_supported
    warning('part1_nozzle_geometry:ExitFlow', ...
        ['Terminal states do not support the assumed uniform exit ', ...
         'region. That region will not be filled automatically.']);
end

%% Interpolate the upper-half Mach field

sample_points = [mesh_x(:), mesh_y(:)];

% Average values at any duplicate sample locations.
[unique_points, ~, point_groups] = unique(sample_points, 'rows');
unique_M = accumarray(point_groups, mesh_M(:), [], @mean);

% Do not extrapolate beyond the characteristic network.
F_M = scatteredInterpolant( ...
    unique_points(:, 1), ...
    unique_points(:, 2), ...
    unique_M, ...
    'natural', ...
    'none');

x_query = linspace(0, x_wall_full(end), plot_nx);
y_query = linspace(0, max(y_wall_full), plot_ny);

[Xq, Yq] = meshgrid(x_query, y_query);

[Mq, inside_nozzle, uniform_mask] = evaluate_mach_field( ...
    F_M, Xq, Yq, ...
    x_wall_full, y_wall_full, ...
    uniform_x, uniform_y, ...
    M_exit, uniform_flow_supported);

%% Reflect the upper-half field across the centerline

% Exclude the first upper-half row from the reflection so that the
% centerline, y = 0, appears only once.
%
% Mach number is a scalar, so its reflected value is unchanged.

y_query_full = [-fliplr(y_query(2:end)), y_query];

[Xq_full, Yq_full] = meshgrid(x_query, y_query_full);

Mq_full = [flipud(Mq(2:end, :)); Mq];

inside_nozzle_full = [ ...
    flipud(inside_nozzle(2:end, :)); ...
    inside_nozzle];

uniform_mask_full = [ ...
    flipud(uniform_mask(2:end, :)); ...
    uniform_mask];

%% Store plotting data

% Preserve the upper-half arrays.
result.Xq = Xq;
result.Yq = Yq;
result.Mq = Mq;

result.inside_nozzle = inside_nozzle;
result.uniform_mask = uniform_mask;
result.uniform_flow_supported = uniform_flow_supported;

% Add the full symmetric contour arrays.
result.Xq_full = Xq_full;
result.Yq_full = Yq_full;
result.Mq_full = Mq_full;

result.inside_nozzle_full = inside_nozzle_full;
result.uniform_mask_full = uniform_mask_full;

result.x_wall_full = x_wall_full;
result.y_wall_full = y_wall_full;

%% Shared plot settings

contour_levels = linspace(1, M_exit, n_colors + 1);

nozzle_length = x_wall_full(end);
maximum_half_height = max(y_wall_full);

% Small margins keep the outline from being clipped.
x_padding = 0.015 * nozzle_length;
y_padding = 0.06 * maximum_half_height;

%% Figure 1: Upper-half geometry and characteristics

figure( ...
    'Name', 'Nozzle Geometry and Characteristics', ...
    'NumberTitle', 'off', ...
    'Color', 'w', ...
    'Position', [100, 100, 1200, 460]);

ax_geometry = axes;
hold(ax_geometry, 'on');

cminus_color = [0.12, 0.38, 0.75];
cplus_color  = [0.65, 0.20, 0.55];

% Representative handles for the legend.
h_cminus = plot(ax_geometry, NaN, NaN, ...
    '-', ...
    'Color', cminus_color, ...
    'LineWidth', 0.7);

h_cplus = plot(ax_geometry, NaN, NaN, ...
    '-', ...
    'Color', cplus_color, ...
    'LineWidth', 0.7);

for k = 1:numel(segments)

    if strcmp(segments(k).type, 'Cminus')
        segment_color = cminus_color;
    else
        segment_color = cplus_color;
    end

    plot(ax_geometry, ...
        segments(k).x, segments(k).y, ...
        '-', ...
        'Color', segment_color, ...
        'LineWidth', 0.7, ...
        'HandleVisibility', 'off');
end

h_wall = plot(ax_geometry, ...
    x_wall_full, y_wall_full, ...
    '-', ...
    'Color', wall_color, ...
    'LineWidth', wall_line_width);

h_centerline = plot(ax_geometry, ...
    [0, nozzle_length], [0, 0], ...
    '--', ...
    'Color', [0.45, 0.45, 0.45], ...
    'LineWidth', 1);

plot(ax_geometry, mesh_x, mesh_y, ...
    '.', ...
    'Color', wall_color, ...
    'MarkerSize', 3, ...
    'HandleVisibility', 'off');

xlabel(ax_geometry, 'x / h_t');
ylabel(ax_geometry, 'y / h_t');

title(ax_geometry, sprintf( ...
    'Mach %.1f Nozzle: Upper-Half Characteristic Network', M_exit));

legend(ax_geometry, ...
    [h_wall, h_cminus, h_cplus, h_centerline], ...
    {'Nozzle wall', 'C^- characteristics', ...
     'C^+ characteristics', 'Centerline'}, ...
    'Location', 'best');

axis(ax_geometry, 'equal');

xlim(ax_geometry, [0, nozzle_length]);
ylim(ax_geometry, [0, 1.05 * maximum_half_height]);

grid(ax_geometry, 'on');
box(ax_geometry, 'on');
set(ax_geometry, 'FontSize', 12);

hold(ax_geometry, 'off');

%% Figure 2: Full symmetric scientific Mach contour

figure( ...
    'Name', 'Full Nozzle Mach Contour', ...
    'NumberTitle', 'off', ...
    'Color', 'w', ...
    'Position', [120, 150, 1200, 520]);

ax_contour = axes;
hold(ax_contour, 'on');

contourf(ax_contour, Xq_full, Yq_full, Mq_full, ...
    contour_levels, ...
    'LineColor', 'none');

colormap(ax_contour, turbo(n_colors));
caxis(ax_contour, [1, M_exit]);

% Draw the full nozzle perimeter.
draw_nozzle_outline( ...
    ax_contour, x_wall_full, y_wall_full, ...
    wall_color, wall_line_width, end_line_width);

% Show the symmetry plane on the scientific plot only.
plot(ax_contour, ...
    [0, nozzle_length], [0, 0], ...
    '--', ...
    'Color', [0.35, 0.35, 0.35], ...
    'LineWidth', 0.8);

cb_contour = colorbar(ax_contour, 'eastoutside');
cb_contour.Label.String = 'Mach number';
cb_contour.Label.FontSize = 13;
cb_contour.FontSize = 11;
cb_contour.Ticks = 1:M_exit;

xlabel(ax_contour, 'x / h_t');
ylabel(ax_contour, 'y / h_t');

title(ax_contour, sprintf( ...
    'Symmetric Nozzle Mach Distribution — Design Exit Mach %.1f', ...
    M_exit));

axis(ax_contour, 'equal');

xlim(ax_contour, [-x_padding, nozzle_length + x_padding]);

ylim(ax_contour, [ ...
    -maximum_half_height - y_padding, ...
     maximum_half_height + y_padding]);

grid(ax_contour, 'off');
box(ax_contour, 'on');
set(ax_contour, 'FontSize', 12);

hold(ax_contour, 'off');

%% Figure 3: Clean full contour with outline and color legend
% No axes, title, grid, or centerline.
% This is a display-only figure. Nothing is saved automatically.

fig_clean = figure( ...
    'Name', 'Full Nozzle Mach Contour — Portfolio View', ...
    'NumberTitle', 'off', ...
    'Color', 'w', ...
    'Position', [150, 200, 1400, 520]);

ax_clean = axes('Parent', fig_clean);
hold(ax_clean, 'on');

contourf(ax_clean, Xq_full, Yq_full, Mq_full, ...
    contour_levels, ...
    'LineColor', 'none');

colormap(ax_clean, turbo(n_colors));
caxis(ax_clean, [1, M_exit]);

draw_nozzle_outline( ...
    ax_clean, x_wall_full, y_wall_full, ...
    wall_color, wall_line_width, end_line_width);

axis(ax_clean, 'equal');

xlim(ax_clean, [-x_padding, nozzle_length + x_padding]);

ylim(ax_clean, [ ...
    -maximum_half_height - y_padding, ...
     maximum_half_height + y_padding]);

axis(ax_clean, 'off');
set(ax_clean, 'Color', 'none');

title(ax_clean, 'Mach 6 Nozzle — Mach Number Distribution', ...
    'FontSize', 18, ...
    'FontWeight', 'bold', ...
    'Color', wall_color, ...
    'Visible', 'on');

cb_clean = colorbar(ax_clean, 'eastoutside');
cb_clean.Label.String = 'Mach number';
cb_clean.Label.FontSize = 14;
cb_clean.FontSize = 12;
cb_clean.Color = [0.15, 0.20, 0.27];
cb_clean.Ticks = 1:M_exit;
cb_clean.Box = 'off';

hold(ax_clean, 'off');

drawnow;

end

%% Local function: Draw the full nozzle outline

function draw_nozzle_outline( ...
    ax, x_wall, y_wall, outline_color, wall_width, end_width)

% Upper physical wall.
plot(ax, x_wall, y_wall, ...
    '-', ...
    'Color', outline_color, ...
    'LineWidth', wall_width, ...
    'HandleVisibility', 'off');

% Lower physical wall, reflected across the centerline.
plot(ax, x_wall, -y_wall, ...
    '-', ...
    'Color', outline_color, ...
    'LineWidth', wall_width, ...
    'HandleVisibility', 'off');

% Inlet and exit lines close the displayed flow-domain outline.
% These are open flow boundaries, not physical solid walls.
plot(ax, ...
    [x_wall(1), x_wall(1)], ...
    [-y_wall(1), y_wall(1)], ...
    '-', ...
    'Color', outline_color, ...
    'LineWidth', end_width, ...
    'HandleVisibility', 'off');

plot(ax, ...
    [x_wall(end), x_wall(end)], ...
    [-y_wall(end), y_wall(end)], ...
    '-', ...
    'Color', outline_color, ...
    'LineWidth', end_width, ...
    'HandleVisibility', 'off');

end

%% Local function: Method-of-characteristics solver

function result = run_moc_nozzle(N, gamma, M_exit, ht, Rt)

nu_exit = PM(M_exit, gamma);
theta_max = nu_exit / 2;

%% Initial expansion points

theta0 = linspace(theta_max / N, theta_max, N);
nu0 = theta0;

M0 = arrayfun(@(nu) inversePM(nu, gamma), nu0);
mu0 = asind(1 ./ M0);

x0 = Rt .* sind(theta0);
y0 = ht + Rt .* (1 - cosd(theta0));

active.x = x0;
active.y = y0;
active.theta = theta0;
active.nu = nu0;
active.M = M0;
active.mu = mu0;

active.Kplus = theta0 - nu0;
active.Kminus = theta0 + nu0;

%% Storage

mesh_x = x0;
mesh_y = y0;
mesh_M = M0;

% The last initial-expansion point starts the designed wall.
previous_wall = get_point(active, N);

x_wall = previous_wall.x;
y_wall = previous_wall.y;
theta_wall = previous_wall.theta;
nu_wall = previous_wall.nu;
M_wall = previous_wall.M;
mu_wall = previous_wall.mu;

x_center = zeros(1, N);
M_center = zeros(1, N);
nu_center = zeros(1, N);

segments = struct('x', {}, 'y', {}, 'type', {});

%% March through the characteristic network

for row = 1:N

    n_active = numel(active.x);

    new = struct( ...
        'x', [], ...
        'y', [], ...
        'theta', [], ...
        'nu', [], ...
        'M', [], ...
        'mu', [], ...
        'Kplus', [], ...
        'Kminus', []);

    %% Centerline point

    A = get_point(active, 1);

    C.theta = 0;
    C.Kminus = A.Kminus;
    C.nu = C.Kminus;
    C.Kplus = C.theta - C.nu;

    C.M = inversePM(C.nu, gamma);
    C.mu = asind(1 / C.M);

    slope_minus = tand( ...
        0.5 * ((A.theta - A.mu) + (C.theta - C.mu)));

    C.x = A.x - A.y / slope_minus;
    C.y = 0;

    new = append_point(new, C);

    x_center(row) = C.x;
    M_center(row) = C.M;
    nu_center(row) = C.nu;

    segments(end + 1) = struct( ...
        'x', [A.x, C.x], ...
        'y', [A.y, C.y], ...
        'type', 'Cminus');

    %% Interior intersections

    for j = 2:n_active

        left = get_point(new, j - 1);
        upper = get_point(active, j);

        P.Kplus = left.Kplus;
        P.Kminus = upper.Kminus;

        P.theta = 0.5 * (P.Kminus + P.Kplus);
        P.nu = 0.5 * (P.Kminus - P.Kplus);

        P.M = inversePM(P.nu, gamma);
        P.mu = asind(1 / P.M);

        slope_plus = tand( ...
            0.5 * ((left.theta + left.mu) + ...
                   (P.theta + P.mu)));

        slope_minus = tand( ...
            0.5 * ((upper.theta - upper.mu) + ...
                   (P.theta - P.mu)));

        P.x = ( ...
            upper.y - left.y + ...
            slope_plus * left.x - ...
            slope_minus * upper.x) / ...
            (slope_plus - slope_minus);

        P.y = left.y + slope_plus * (P.x - left.x);

        new = append_point(new, P);

        segments(end + 1) = struct( ...
            'x', [left.x, P.x], ...
            'y', [left.y, P.y], ...
            'type', 'Cplus');

        segments(end + 1) = struct( ...
            'x', [upper.x, P.x], ...
            'y', [upper.y, P.y], ...
            'type', 'Cminus');
    end

    %% Wall point

    last_interior = get_point(new, numel(new.x));

    W.Kplus = last_interior.Kplus;
    W.Kminus = nu_exit;

    W.theta = 0.5 * (W.Kminus + W.Kplus);
    W.nu = 0.5 * (W.Kminus - W.Kplus);

    W.M = inversePM(W.nu, gamma);
    W.mu = asind(1 / W.M);

    slope_plus = tand( ...
        0.5 * ((last_interior.theta + last_interior.mu) + ...
               (W.theta + W.mu)));

    slope_wall = tand( ...
        0.5 * (previous_wall.theta + W.theta));

    W.x = ( ...
        previous_wall.y - last_interior.y + ...
        slope_plus * last_interior.x - ...
        slope_wall * previous_wall.x) / ...
        (slope_plus - slope_wall);

    W.y = previous_wall.y + ...
        slope_wall * (W.x - previous_wall.x);

    segments(end + 1) = struct( ...
        'x', [last_interior.x, W.x], ...
        'y', [last_interior.y, W.y], ...
        'type', 'Cplus');

    x_wall(end + 1) = W.x;
    y_wall(end + 1) = W.y;
    theta_wall(end + 1) = W.theta;
    nu_wall(end + 1) = W.nu;
    M_wall(end + 1) = W.M;
    mu_wall(end + 1) = W.mu;

    previous_wall = W;

    %% Store computed points

    mesh_x = [mesh_x, new.x, W.x];
    mesh_y = [mesh_y, new.y, W.y];
    mesh_M = [mesh_M, new.M, W.M];

    %% Prepare the next row

    fields = fieldnames(active);

    for k = 1:numel(fields)
        field = fields{k};
        row_values = new.(field);
        active.(field) = [row_values(2:end), W.(field)];
    end
end

%% Return results

result.x0 = x0;
result.y0 = y0;
result.theta0 = theta0;
result.nu0 = nu0;
result.M0 = M0;
result.mu0 = mu0;

result.x_wall = x_wall;
result.y_wall = y_wall;
result.theta_wall = theta_wall;
result.nu_wall = nu_wall;
result.M_wall = M_wall;
result.mu_wall = mu_wall;

result.x_center = x_center;
result.M_center = M_center;
result.nu_center = nu_center;

result.mesh_x = mesh_x;
result.mesh_y = mesh_y;
result.mesh_M = mesh_M;

result.segments = segments;

result.nu_exit = nu_exit;
result.theta_max = theta_max;

% Retain the original half-height fields for compatibility.
result.exit_height = y_wall(end);
result.nozzle_height = y_wall(end);
result.nozzle_length = x_wall(end);

% Full height is provided separately.
result.full_exit_height = 2 * y_wall(end);

% For a planar nozzle of constant width, the area ratio equals the
% height ratio. Using half-heights gives the same ratio as full heights.
result.area_ratio = y_wall(end) / ht;

result.area_ratio_exact = (1 / M_exit) * ...
    ((2 / (gamma + 1)) * ...
    (1 + ((gamma - 1) / 2) * M_exit^2)) ^ ...
    ((gamma + 1) / (2 * (gamma - 1)));

end

%% Local function: Evaluate the upper-half Mach field

function [Mq, inside_nozzle, uniform_mask] = evaluate_mach_field( ...
    F_M, Xq, Yq, ...
    x_wall, y_wall, ...
    uniform_x, uniform_y, ...
    M_exit, fill_uniform)

% Interpolate the upper-wall height at each query location.
wall_height = interp1( ...
    x_wall, y_wall, Xq, 'linear', 'extrap');

inside_nozzle = ...
    Xq >= x_wall(1) & ...
    Xq <= x_wall(end) & ...
    Yq >= 0 & ...
    Yq <= wall_height;

% Interpolate within the characteristic network.
Mq = F_M(Xq, Yq);

% Fill only the identified uniform exit-flow triangle.
uniform_mask = false(size(Xq));

if fill_uniform
    uniform_mask = inside_nozzle & ...
        inpolygon(Xq, Yq, uniform_x, uniform_y);

    Mq(uniform_mask) = M_exit;
end

% Hide all points outside the nozzle.
Mq(~inside_nozzle) = NaN;

end

%% Local function: Prandtl-Meyer angle in degrees

function nu = PM(M, gamma)

nu = sqrt((gamma + 1) / (gamma - 1)) .* ...
    atand(sqrt(((gamma - 1) / (gamma + 1)) .* (M.^2 - 1))) ...
    - atand(sqrt(M.^2 - 1));

end

%% Local function: Invert the Prandtl-Meyer relation

function M = inversePM(nu, gamma)

if abs(nu) < 1e-12
    M = 1;
    return;
end

M = fzero(@(mach) PM(mach, gamma) - nu, [1, 100]);

end

%% Local function: Append a point to a characteristic row

function data = append_point(data, point)

fields = fieldnames(data);

for k = 1:numel(fields)
    field = fields{k};
    data.(field)(end + 1) = point.(field);
end

end

%% Local function: Extract one point from a characteristic row

function point = get_point(data, index)

fields = fieldnames(data);

for k = 1:numel(fields)
    field = fields{k};
    point.(field) = data.(field)(index);
end

end