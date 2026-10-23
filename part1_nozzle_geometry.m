function result = part1_nozzle_geometry(N, make_plots)
% RECONSTRUCTED COMPONENT: Geometry, characteristics, and Mach map.
% Extracted from the completed integrated project.
%
% Normal execution:
%   result = part1_nozzle_geometry;
%
% Numerical solution only, used by Parts 2 and 3:
%   result = part1_nozzle_geometry(50, false);

if nargin < 1
    N = 50;
end

if nargin < 2
    make_plots = true;
end

%% DESIGN INPUTS

gamma = 1.4;
M_exit = 6;

ht = 1;
Rt = ht;

%% CALCULATE NOZZLE SOLUTION

result = run_moc_nozzle(N, gamma, M_exit, ht, Rt);

% Pass design inputs to the other project components.
result.gamma = gamma;
result.M_exit = M_exit;
result.ht = ht;
result.Rt = Rt;
result.N = N;

% Parts 2 and 3 request the numerical solution without plots.
if ~make_plots
    return
end

%% EXTRACT RESULTS

x0 = result.x0;
y0 = result.y0;

x_wall = result.x_wall;
y_wall = result.y_wall;
theta_wall = result.theta_wall;
M_wall = result.M_wall;

x_center = result.x_center;
M_center = result.M_center;

mesh_x = result.mesh_x;
mesh_y = result.mesh_y;
mesh_M = result.mesh_M;

seg_x1 = result.seg_x1;
seg_y1 = result.seg_y1;
seg_x2 = result.seg_x2;
seg_y2 = result.seg_y2;
seg_type = result.seg_type;

Ae_At_expected = (1/M_exit) * ...
    ((2/(gamma+1)) * ...
    (1 + ((gamma-1)/2)*M_exit^2)) ...
    ^((gamma+1)/(2*(gamma-1)));

L_triangle = x_wall(end) - x_center(end);

%% PRINT RESULTS

fprintf('\n================ NOZZLE GEOMETRY ================\n')
fprintf('Number of characteristics:              %d\n', N)
fprintf('Final centerline Mach number:           %.6f\n', M_center(end))
fprintf('Final wall Mach number:                 %.6f\n', M_wall(end))
fprintf('Target exit Mach number:                %.6f\n', M_exit)
fprintf('Final wall angle:                       %.8f deg\n', theta_wall(end))
fprintf('Exit height ratio, y_e/h_t:             %.6f\n', y_wall(end)/ht)
fprintf('Expected area ratio, A_e/A_t:           %.6f\n', Ae_At_expected)
fprintf('Height ratio error:                     %.6f\n', ...
    y_wall(end)/ht - Ae_At_expected)
fprintf('Nozzle length, x_e/h_t:                 %.6f\n', x_wall(end))
fprintf('Uniform-flow triangle length, L/h_t:    %.6f\n\n', L_triangle)

%% FULL WALL CONTOUR INCLUDING THROAT ARC

x_wall_full = [0, x0, x_wall(2:end)];
y_wall_full = [ht, y0, y_wall(2:end)];

[x_wall_unique, unique_idx] = unique(x_wall_full, 'stable');
y_wall_unique = y_wall_full(unique_idx);

%% PLOT 1: NOZZLE CONTOUR AND CHARACTERISTICS

figure
hold on
grid on
axis equal

h_wall = plot(x_wall_full, y_wall_full, ...
    'r-', 'LineWidth', 2.5);

h_arc = plot(x0, y0, ...
    'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

h_cminus = [];
h_cplus = [];

for s = 1:length(seg_type)

    if seg_type(s) == "Cminus"

        h = plot([seg_x1(s), seg_x2(s)], ...
            [seg_y1(s), seg_y2(s)], ...
            'b', 'LineWidth', 0.6);

        if isempty(h_cminus)
            h_cminus = h;
        end

    elseif seg_type(s) == "Cplus"

        h = plot([seg_x1(s), seg_x2(s)], ...
            [seg_y1(s), seg_y2(s)], ...
            'm', 'LineWidth', 0.6);

        if isempty(h_cplus)
            h_cplus = h;
        end
    end
end

h_center = plot(x_center, zeros(size(x_center)), ...
    'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

yline(0, 'k--', 'Centerline')

xlabel('x / h_t')
ylabel('y / h_t')
title(sprintf('Nozzle Expansion Contour and Characteristics, N = %d', N))

legend([h_wall h_arc h_cminus h_cplus h_center], ...
    {'Nozzle wall contour', ...
     'Initial throat points', ...
     'C^- characteristic', ...
     'C^+ characteristic', ...
     'Centerline points'}, ...
    'Location', 'best')

%% INTERPOLATION GRID FOR MACH MAP

nx = 350;
ny = 200;

xq = linspace(0, x_wall_full(end), nx);
yq = linspace(0, max(y_wall_full), ny);

[Xq, Yq] = meshgrid(xq, yq);

y_wall_interp = interp1( ...
    x_wall_unique, y_wall_unique, xq, 'linear', 'extrap');

Ywallq = repmat(y_wall_interp, ny, 1);

inside_nozzle = Yq >= 0 & Yq <= Ywallq;

%% PLOT 2: MACH CONTOUR MAP

F_M = scatteredInterpolant( ...
    mesh_x(:), mesh_y(:), mesh_M(:), 'natural', 'none');

Mq = F_M(Xq, Yq);
Mq(~inside_nozzle) = NaN;

figure
hold on

contourf(Xq, Yq, Mq, 40, 'LineColor', 'none')
colorbar
colormap turbo

plot(x_wall_full, y_wall_full, 'k-', 'LineWidth', 2)
plot([0 x_wall_full(end)], [0 0], 'k--', 'LineWidth', 1)

axis equal
grid on

xlabel('x / h_t')
ylabel('y / h_t')
title(sprintf('Mach Number Contour Map, N = %d', N))
caxis([1 M_exit])

end

%% MAIN SOLVER FUNCTION

function result = run_moc_nozzle(N, gamma, M_exit, ht, Rt)

nu_exit = PM(M_exit, gamma);
theta_max = nu_exit/2;

theta0 = linspace(theta_max/N, theta_max, N);
nu0 = theta0;

x0 = Rt*sind(theta0);
y0 = ht + Rt*(1 - cosd(theta0));

M0 = zeros(1, N);
mu0 = zeros(1, N);

for i = 1:N
    M0(i) = inversePM(nu0(i), gamma);
    mu0(i) = asind(1/M0(i));
end

Kplus0 = theta0 - nu0;
Kminus0 = theta0 + nu0;

active.x = x0;
active.y = y0;
active.theta = theta0;
active.nu = nu0;
active.M = M0;
active.mu = mu0;
active.Kplus = Kplus0;
active.Kminus = Kminus0;

seg_x1 = [];
seg_y1 = [];
seg_x2 = [];
seg_y2 = [];
seg_type = strings(0);

mesh_x = x0;
mesh_y = y0;
mesh_M = M0;

x_center = [];
M_center = [];
nu_center = [];

x_wall = x0(end);
y_wall = y0(end);
theta_wall = theta0(end);
nu_wall = nu0(end);
M_wall = M0(end);

for row = 1:N

    nActive = length(active.x);

    new.x = [];
    new.y = [];
    new.theta = [];
    new.nu = [];
    new.M = [];
    new.mu = [];
    new.Kplus = [];
    new.Kminus = [];

    %% 1. Centerline point from first active C-minus

    A = get_point(active, 1);

    C.theta = 0;
    C.nu = A.Kminus;
    C.M = inversePM(C.nu, gamma);
    C.mu = asind(1/C.M);

    C.Kplus = C.theta - C.nu;
    C.Kminus = C.theta + C.nu;

    m_minus = tand(0.5 * ...
        ((A.theta - A.mu) + (C.theta - C.mu)));

    C.x = A.x - A.y/m_minus;
    C.y = 0;

    new = append_point(new, C);

    x_center(end+1) = C.x;
    M_center(end+1) = C.M;
    nu_center(end+1) = C.nu;

    mesh_x(end+1) = C.x;
    mesh_y(end+1) = C.y;
    mesh_M(end+1) = C.M;

    seg_x1(end+1) = A.x;
    seg_y1(end+1) = A.y;
    seg_x2(end+1) = C.x;
    seg_y2(end+1) = C.y;
    seg_type(end+1) = "Cminus";

    %% 2. Interior points

    for j = 2:nActive

        left = get_point(new, j-1);
        upper = get_point(active, j);

        P.Kplus = left.Kplus;
        P.Kminus = upper.Kminus;

        P.theta = 0.5*(P.Kminus + P.Kplus);
        P.nu = 0.5*(P.Kminus - P.Kplus);

        P.M = inversePM(P.nu, gamma);
        P.mu = asind(1/P.M);

        m_plus = tand(0.5 * ...
            ((left.theta + left.mu) + (P.theta + P.mu)));

        m_minus = tand(0.5 * ...
            ((upper.theta - upper.mu) + (P.theta - P.mu)));

        P.x = (upper.y - left.y ...
            + m_plus*left.x - m_minus*upper.x) ...
            / (m_plus - m_minus);

        P.y = left.y + m_plus*(P.x - left.x);

        new = append_point(new, P);

        mesh_x(end+1) = P.x;
        mesh_y(end+1) = P.y;
        mesh_M(end+1) = P.M;

        seg_x1(end+1) = left.x;
        seg_y1(end+1) = left.y;
        seg_x2(end+1) = P.x;
        seg_y2(end+1) = P.y;
        seg_type(end+1) = "Cplus";

        seg_x1(end+1) = upper.x;
        seg_y1(end+1) = upper.y;
        seg_x2(end+1) = P.x;
        seg_y2(end+1) = P.y;
        seg_type(end+1) = "Cminus";
    end

    %% 3. Wall point from tangency and compatibility

    lastInterior = get_point(new, length(new.x));

    W.Kplus = lastInterior.Kplus;
    W.Kminus = nu_exit;

    W.theta = 0.5*(W.Kminus + W.Kplus);
    W.nu = 0.5*(W.Kminus - W.Kplus);

    W.M = inversePM(W.nu, gamma);
    W.mu = asind(1/W.M);

    prevWall.x = x_wall(end);
    prevWall.y = y_wall(end);
    prevWall.theta = theta_wall(end);

    m_wall = tand(0.5*(prevWall.theta + W.theta));

    m_plus = tand(0.5 * ...
        ((lastInterior.theta + lastInterior.mu) ...
        + (W.theta + W.mu)));

    W.x = (lastInterior.y - prevWall.y ...
        + m_wall*prevWall.x - m_plus*lastInterior.x) ...
        / (m_wall - m_plus);

    W.y = prevWall.y + m_wall*(W.x - prevWall.x);

    x_wall(end+1) = W.x;
    y_wall(end+1) = W.y;
    theta_wall(end+1) = W.theta;
    nu_wall(end+1) = W.nu;
    M_wall(end+1) = W.M;

    mesh_x(end+1) = W.x;
    mesh_y(end+1) = W.y;
    mesh_M(end+1) = W.M;

    seg_x1(end+1) = lastInterior.x;
    seg_y1(end+1) = lastInterior.y;
    seg_x2(end+1) = W.x;
    seg_y2(end+1) = W.y;
    seg_type(end+1) = "Cplus";

    %% 4. Advance active row

    if row < N

        active.x = [new.x(2:end), W.x];
        active.y = [new.y(2:end), W.y];
        active.theta = [new.theta(2:end), W.theta];
        active.nu = [new.nu(2:end), W.nu];
        active.M = [new.M(2:end), W.M];
        active.mu = [new.mu(2:end), W.mu];
        active.Kplus = [new.Kplus(2:end), W.Kplus];
        active.Kminus = [new.Kminus(2:end), W.Kminus];
    end
end

result.x0 = x0;
result.y0 = y0;

result.x_wall = x_wall;
result.y_wall = y_wall;
result.theta_wall = theta_wall;
result.nu_wall = nu_wall;
result.M_wall = M_wall;

result.x_center = x_center;
result.M_center = M_center;
result.nu_center = nu_center;

result.mesh_x = mesh_x;
result.mesh_y = mesh_y;
result.mesh_M = mesh_M;

result.seg_x1 = seg_x1;
result.seg_y1 = seg_y1;
result.seg_x2 = seg_x2;
result.seg_y2 = seg_y2;
result.seg_type = seg_type;

end

%% HELPER FUNCTIONS

function nu = PM(M, gamma)

if nargin < 2
    gamma = 1.4;
end

if any(M < 1)
    error('Mach number must be >= 1.')
end

term1 = sqrt((gamma+1)/(gamma-1));

term2 = atan(sqrt( ...
    ((gamma-1)/(gamma+1)).*(M.^2 - 1)));

term3 = atan(sqrt(M.^2 - 1));

nu = rad2deg(term1.*term2 - term3);

end

function M = inversePM(nu, gamma)

if nargin < 2
    gamma = 1.4;
end

if any(nu < -1e-10)
    error('Prandtl-Meyer angle must be nonnegative.')
end

if abs(nu) < 1e-10
    M = 1;
    return
end

f = @(Mach) PM(Mach, gamma) - nu;
M = fzero(f, [1.000001, 100]);

end

function data = append_point(data, P)

data.x(end+1) = P.x;
data.y(end+1) = P.y;
data.theta(end+1) = P.theta;
data.nu(end+1) = P.nu;
data.M(end+1) = P.M;
data.mu(end+1) = P.mu;

data.Kplus(end+1) = P.theta - P.nu;
data.Kminus(end+1) = P.theta + P.nu;

end

function P = get_point(data, index)

P.x = data.x(index);
P.y = data.y(index);
P.theta = data.theta(index);
P.nu = data.nu(index);
P.M = data.M(index);
P.mu = data.mu(index);
P.Kplus = data.Kplus(index);
P.Kminus = data.Kminus(index);

end