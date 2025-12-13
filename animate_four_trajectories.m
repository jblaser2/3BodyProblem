% animate_four_trajectories.m
% Integrate four CR3BP trajectories and animate them in the inertial frame
% so the two primaries (m1,m2) appear to move on circular orbits.
%
% Requires: cr3bp_ode.m, lagrange_points.m, perturb_hook.m (optional), jacobi_const.m
% Save file in same folder and run: >> animate_four_trajectories

clear; close all; clc;

%% PARAMETERS
mu = 0.0121505856;           % Earth-Moon like
perturb_eps = 0.0;           % small perturbation amplitude (0 => no perturb)
tspan = [0, 200];            % nondimensional time interval
solver_opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

% time grid for synchronized animation (dense enough)
t_anim = linspace(tspan(1), tspan(2), 3500)';   % column

%% FOUR initial conditions (in rotating frame)
% 1: near L1 (slightly to left)
[Lx, Ly] = lagrange_points(mu);
ic1 = [Lx(1) - 0.02; 0; 0; 0.03];

% 2: near L2 (slightly to right)
ic2 = [Lx(2) + 0.02; 0; 0; -0.02];

% 3: near L4 (triangle)
ic3 = [Lx(4) + 0.02; Ly(4) - 0.02; 0; 0];

% 4: a passing trajectory between primaries (transfer-like)
ic4 = [0.2; 0.05; -0.02; 0.02];

initial_conds = [ic1, ic2, ic3, ic4];

numTraj = size(initial_conds,2);
solT = cell(1,numTraj);
solY = cell(1,numTraj);

%% INTEGRATE each trajectory separately (adaptive times)
fprintf('Integrating %d trajectories with ode45...\n', numTraj);
for k = 1:numTraj
    y0 = initial_conds(:,k);
    [t,Y] = ode45(@(tt,yy) cr3bp_ode(tt,yy,mu,@perturb_hook,perturb_eps), tspan, y0, solver_opts);
    solT{k} = t(:);
    solY{k} = Y;  % columns: x, y, xd, yd
    fprintf('  Traj %d integrated: %d steps, t_end=%.3f\n', k, length(t), t(end));
end

%% Interpolate each solution onto common animation time grid
Y_interp = zeros(length(t_anim), 4, numTraj); % time x state x traj
for k=1:numTraj
    % interp1 requires monotonically increasing solT{k}
    Y_interp(:,:,k) = interp1(solT{k}, solY{k}, t_anim, 'linear', 'extrap');
end

%% Convert rotating-frame positions (x,y) to inertial-frame (X,Y)
% Rotation angle = theta(t) = t (since nondimensional omega = 1)
theta = t_anim;        % angular rotation of the rotating frame
cosT = cos(theta); sinT = sin(theta);

% primaries in rotating frame are at (-mu,0) and (1-mu,0).
prim1_rot = [-mu; 0];
prim2_rot = [1-mu; 0];

% precompute inertial positions for primaries
P1_inertial = [ prim1_rot(1)*cosT - prim1_rot(2)*sinT,  ...
                prim1_rot(1)*sinT + prim1_rot(2)*cosT ]; % Nx2
P2_inertial = [ prim2_rot(1)*cosT - prim2_rot(2)*sinT,  ...
                prim2_rot(1)*sinT + prim2_rot(2)*cosT ];

% trajectories inertial coordinates
Traj_inertial = zeros(length(t_anim),2,numTraj);
for k=1:numTraj
    xr = Y_interp(:,1,k);  yr = Y_interp(:,2,k);
    % rotation: X = x*cos(t) - y*sin(t); Y = x*sin(t) + y*cos(t)
    Traj_inertial(:,1,k) = xr .* cosT - yr .* sinT;
    Traj_inertial(:,2,k) = xr .* sinT + yr .* cosT;
end

%% Prepare figure for animation
figure('Name','Four CR3BP trajectories (inertial frame)','NumberTitle','off','Color','w',...
       'Position',[100 80 1000 700]);
ax = axes('XLim',[-1.5 1.5],'YLim',[-1.2 1.2],'NextPlot','add');
axis equal; grid on;
xlabel('X (inertial)'); ylabel('Y (inertial)');
title('Four CR3BP trajectories with moving primaries (inertial frame)');

% plot circular orbit guides of primaries for context
r1 = norm(prim1_rot);
r2 = norm(prim2_rot);
thc = linspace(0,2*pi,400);
plot(r1*cos(thc), r1*sin(thc), ':','LineWidth',0.8);
plot(r2*cos(thc), r2*sin(thc), ':','LineWidth',0.8);

% create objects for trails and moving markers
trailLen = 300;  % max number of points in visible trailing tail
colors = lines(numTraj);
hTrail = gobjects(numTraj,1);
hPoint = gobjects(numTraj,1);
for k=1:numTraj
    hTrail(k) = plot(ax, NaN, NaN, '-', 'LineWidth', 1.2, 'Color', colors(k,:));
    hPoint(k) = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', colors(k,:), 'MarkerSize',6, ...
                     'MarkerEdgeColor','k');
end
% primaries markers
hP1 = plot(ax, NaN, NaN, 'ro', 'MarkerSize',10, 'MarkerFaceColor','r');
hP2 = plot(ax, NaN, NaN, 'ko', 'MarkerSize',8, 'MarkerFaceColor','k');

legend_entries = arrayfun(@(k) sprintf('traj %d',k), 1:numTraj, 'UniformOutput', false);
legend([hPoint; hP1; hP2], [legend_entries, {'Primary 1','Primary 2'}], 'Location','northeastoutside');

%% Animation loop
fprintf('Animating...\n');
drawnow;
for i = 1:length(t_anim)
    % update primaries
    set(hP1, 'XData', P1_inertial(i,1), 'YData', P1_inertial(i,2));
    set(hP2, 'XData', P2_inertial(i,1), 'YData', P2_inertial(i,2));
    % update each trajectory: trail + moving point
    for k=1:numTraj
        Xi = Traj_inertial(1:max(1,i),1,k);
        Yi = Traj_inertial(1:max(1,i),2,k);
        % show only trailing part
        startIdx = max(1, i-trailLen+1);
        set(hTrail(k), 'XData', Traj_inertial(startIdx:i,1,k), 'YData', Traj_inertial(startIdx:i,2,k));
        set(hPoint(k), 'XData', Traj_inertial(i,1,k), 'YData', Traj_inertial(i,2,k));
    end
    % optional: update title with time
    tt = t_anim(i);
    title(sprintf('Four CR3BP trajectories (inertial frame), t = %.2f', tt));
    drawnow limitrate;
end
fprintf('Animation finished.\n');
