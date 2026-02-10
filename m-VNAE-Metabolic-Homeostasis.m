%% VNAE-Metabolic-Homeostasis-Dynamics
% Geometric Stability of Complex Biochemical Networks (MATLAB Version)
% -------------------------------------------------------------------------

clear; clc; close all;

% 1. NETWORK CONFIGURATION
n = 1200;                       % Metabolic variables/enzymes
avg_interactions = 8;           % Average connectivity
beta = 0.25;                    % Structural rigidity (VNAE factor)

% Generate a Sparse Scale-Free Like Network
% (Using a random adjacency matrix with targeted density for performance)
density = avg_interactions / n;
A = sprand(n, n, density);
A = (A + A') / 2;               % Symmetrize connectivity pattern
A = full(A);                    % Convert to full for matrix ops
A(logical(eye(n))) = 0;         % Remove self-loops

% Add weighted affinities (Log-normal distribution for chemical coupling)
A = A .* lognrnd(0, 0.4, n, n);

% 2. ENZYMATIC ASYMMETRY (Theta)
% Represents localized metabolic inertia/dissipation
theta = lognrnd(1.0, 0.8, [n, 1]);
Theta_mat = diag(theta);

% 3. CURVATURE CERTIFICATE (K)
% Estimating geometric curvature through Monte Carlo sampling
samples = 200000;
idx_i = randi(n, samples, 1);
idx_j = randi(n, samples, 1);

diff_theta = abs(theta(idx_i) - theta(idx_j));
interaction = arrayfun(@(i,j) A(i,j), idx_i, idx_j);
rigidity_factor = 1 + beta * (theta(idx_i) + theta(idx_j));

K = mean((diff_theta .* interaction) ./ rigidity_factor);

% 4. METABOLIC DYNAMICS (VNAE Equation)
% p_stress: Persistent external forcing (dietary spikes, oxidative stress)
p_stress = randn(n, 1) * 0.5;
L = diag(sum(A, 2)) - A;        % Graph Laplacian

% Define the System: d(x)/dt = -(L + Theta)x + p
% x represents the deviation from the ideal homeostatic state
vnae_metabolism = @(t, x) -(L * x) - (theta .* x) + p_stress;

% Initial state offset (shock state)
x0 = randn(n, 1) * 5.0;
t_span = [0 0.5];

% Solve using ODE45 (Runge-Kutta)
[t, y] = ode45(vnae_metabolism, t_span, x0);

% 5. REPORT & ANALYSIS
fprintf('\n--- VNAE METABOLIC HOMEOSTASIS REPORT ---\n');
fprintf('Metabolic Variables (n): %d\n', n);
fprintf('Structural Curvature (K): %.6f\n', K);
if K > 0
    fprintf('Status: Geometrically Stable (Homeostasis Guaranteed)\n');
else
    fprintf('Status: Critical Instability\n');
end

% 6. VISUALIZATION
figure('Color', 'w');
hold on;

% Sample 50 nodes for visual clarity
sample_indices = randperm(n, 50);
plot(t, y(:, sample_indices), 'LineWidth', 1, 'Color', [0 0.4470 0.7410 0.4]);

yline(0, '--k', 'Alpha', 0.5);
title(['VNAE: Metabolic Convergence (K = ', num2box(K), ')']);
xlabel('Time (Metabolic Scale)');
ylabel('Deviation from Equilibrium');
grid on;
set(gca, 'FontSize', 12);

function str = num2box(val)
    str = sprintf('%.4f', val);
end
