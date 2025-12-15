% 3 Body Problem -- Modified

% Main script for the planar Circular Restricted Three-Body Problem (CR3BP)
% - calls modular functions: cr3bp_ode, lagrange_points, jacobi_const, ...
% - integrates, plots trajectory, Jacobi const, and Poincare section

% This script models a system with 2 large masses (Earth and Moon). The
% system is examined from the rotating reference of these 2 bodies to see
% how a 3rd 'massless' particle

clear; close all; clc;

%% PARAMETERS 
mu = 0.0121505856;    % mass ratio m2/(m1+m2) : Earth-Moon ~0.01215
perturb_eps = 0.0;    % perturbation amplitude (set >0 to enable)
tspan = [0, 100];     % nondimensional time (rotating frame units)
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'Events',@yCrossEvent);

%% Lagrange points
[Lx, Ly] = lagrange_points(mu);

%% Define 5 interesting initial conditions (columns)
% Small offsets from each Lagrange point
eps0 = 0.01;

ICs = [ ...
    Lx(1)+eps0,  Lx(2)+eps0,  Lx(3)-eps0,  Lx(4),  Lx(5);   % x
    Ly(1),      Ly(2),      Ly(3),      Ly(4)+eps0,  Ly(5)-eps0;   % y
    0,          0,          0,          0,          0;           % xd
    0,          0,          0,          0,          0            % yd
];

labels = {'Near L1','Near L2','Near L3','Near L4','Near L5'};
colors = lines(5);

%% Storage
Traj = cell(1,5);
Time = cell(1,5);
Jacobi = cell(1,5);

%% Integrate all trajectories
for k = 1:5
    y0_state = ICs(:,k);
    [t, Y] = ode45(@(t,y) cr3bp_ode(t,y,mu,@perturb_hook,perturb_eps), ...
                    tspan, y0_state, opts);

    Time{k} = t;
    Traj{k} = Y;
    Jacobi{k} = jacobi_const(Y(:,1), Y(:,2), Y(:,3), Y(:,4), mu);
end

%% Plot trajectories in rotating frame
figure('Name','CR3BP trajectories (rotating frame)', ...
       'NumberTitle','off','Position',[100 100 1000 420]);

subplot(1,2,1); hold on;
for k = 1:5
    plot(Traj{k}(:,1), Traj{k}(:,2), 'LineWidth',1.4, ...
         'Color', colors(k,:), 'DisplayName', labels{k});
end

% Primaries
plot(-mu, 0, 'ro','MarkerFaceColor','r','MarkerSize',10,'DisplayName','Primary 1');
plot(1-mu, 0, 'ko','MarkerFaceColor','k','MarkerSize',8,'DisplayName','Primary 2');

% Lagrange points
plot(Lx, Ly, 'g+','MarkerSize',10,'LineWidth',1.5,'DisplayName','Lagrange points');

axis equal; grid on;
xlabel('x'); ylabel('y');
title('Particle trajectories (rotating frame)');
legend('Location','bestoutside');

% Plot Jacobi Constants
subplot(1,2,2); hold on;
for k = 1:5
    plot(Time{k}, Jacobi{k}, 'LineWidth',1.4, ...
         'Color', colors(k,:), 'DisplayName', labels{k});
end

xlabel('t'); ylabel('Jacobi constant C');
title('Jacobi constant along each trajectory');
grid on;
legend('Location','best');

%% Phase space plots
figure('Name','Phase space trajectories (CR3BP)', ...
       'NumberTitle','off','Position',[120 120 1000 700]);

% ---- x vs xdot ----
subplot(2,1,1); hold on;
for k = 1:5
    plot(Traj{k}(:,1), Traj{k}(:,3), ...
         'LineWidth',1.3, ...
         'Color', colors(k,:), ...
         'DisplayName', labels{k});
end
xlabel('x'); ylabel('xdot');
title('Phase space: x vs xdot');
grid on;
legend('Location','best');

% ---- y vs ydot ----
subplot(2,1,2); hold on;
for k = 1:5
    plot(Traj{k}(:,2), Traj{k}(:,4), ...
         'LineWidth',1.3, ...
         'Color', colors(k,:), ...
         'DisplayName', labels{k});
end
xlabel('y'); ylabel('ydot');
title('Phase space: y vs ydot');
grid on;
legend('Location','best');


%% Phase space plots: one figure per trajectory
for k = 1:5

    figure('Name', ['Phase Space - ' labels{k}], ...
           'NumberTitle','off', ...
           'Position',[150+40*k 120 700 600]);

    % --- x vs xdot ---
    subplot(2,1,1);
    plot(Traj{k}(:,1), Traj{k}(:,3), 'LineWidth',1.4);
    grid on;
    xlabel('x');
    ylabel('xdot');
    title(['Phase Space: x vs xdot  (' labels{k} ')']);

    % mark initial condition
    hold on;
    plot(Traj{k}(1,1), Traj{k}(1,3), 'ro', 'MarkerFaceColor','r');
    legend('trajectory','initial condition','Location','best');

    % --- y vs ydot ---
    subplot(2,1,2);
    plot(Traj{k}(:,2), Traj{k}(:,4), 'LineWidth',1.4);
    grid on;
    xlabel('y');
    ylabel('ydot');
    title(['Phase Space: y vs ydot  (' labels{k} ')']);

    % mark initial condition
    hold on;
    plot(Traj{k}(1,2), Traj{k}(1,4), 'ro', 'MarkerFaceColor','r');
    legend('trajectory','initial condition','Location','best');

end


%% FAST animated x-y trajectories (rotating frame)
fps = 20;
Nframes = 250;           % <<< key speed control
delay_time = 1 / fps;

for k = 1:5

    X = Traj{k}(:,1);
    Y = Traj{k}(:,2);

    % Downsample indices
    idx = round(linspace(1, length(X), Nframes));

    % Create figure (smaller = faster)
    fig = figure('Color','w','Position',[300 200 500 450]);
    axis equal;
    hold on;
    xlim([-1.5 1.5]);
    ylim([-1.2 1.2]);

    % Static elements (draw once)
    plot(-mu, 0, 'ro','MarkerFaceColor','r','MarkerSize',8);
    plot(1-mu, 0, 'ko','MarkerFaceColor','k','MarkerSize',6);
    plot(Lx, Ly, 'g+','MarkerSize',8,'LineWidth',1.2);

    xlabel('x'); ylabel('y');
    title(['Trajectory: ' labels{k}]);

    % Graphics handles
    hTrail = plot(NaN,NaN,'b-','LineWidth',1.5);
    hDot   = plot(NaN,NaN,'bo','MarkerFaceColor','b','MarkerSize',6);

    gif_name = sprintf('trajectory_fast_%d.gif', k);

    for j = 1:Nframes
        i = idx(j);

        % Update trajectory
        set(hTrail,'XData',X(idx(1:j)),'YData',Y(idx(1:j)));
        set(hDot,  'XData',X(i),      'YData',Y(i));

        drawnow limitrate;

        % Capture frame
        frame = getframe(fig);
        img = frame2im(frame);
        [imind, cm] = rgb2ind(img, 256);

        if j == 1
            imwrite(imind, cm, gif_name, 'gif', ...
                    'Loopcount', inf, 'DelayTime', delay_time);
        else
            imwrite(imind, cm, gif_name, 'gif', ...
                    'WriteMode','append', 'DelayTime', delay_time);
        end
    end

    close(fig);
    fprintf('Saved %s\n', gif_name);
end


%% Poincaré section for ONE trajectory (y = 0, ydot > 0)

k = 1;   % choose trajectory (e.g. near L4)

Ytraj = Traj{k};   % Nx4: [x, y, xdot, ydot]

xcross = [];
xdotcross = [];

for n = 1:size(Ytraj,1)-1
    y1 = Ytraj(n,2);
    y2 = Ytraj(n+1,2);

    % upward crossing through y = 0
    if (y1 <= 0) && (y2 > 0)

        % linear interpolation
        alpha = -y1 / (y2 - y1);

        x_at = Ytraj(n,1) + alpha*(Ytraj(n+1,1) - Ytraj(n,1));
        xd_at = Ytraj(n,3) + alpha*(Ytraj(n+1,3) - Ytraj(n,3));

        xcross(end+1,1) = x_at;
        xdotcross(end+1,1) = xd_at;
    end
end

% Plot Poincaré section
figure('Name',['Poincare Section - ' labels{k}], ...
       'NumberTitle','off','Position',[300 200 600 500]);

plot(xcross, xdotcross, 'k.', 'MarkerSize',8);
grid on;

xlabel('x');
ylabel('xdot');

title(['Poincaré Section (y=0, ydot>0) -- ' labels{k}]);

figure;
plot(Ytraj(:,2)); % Plot y values over time
xlabel('Time step');
ylabel('y');
title('y values over time');
grid on;

%% Choose ONE trajectory for detailed analysis
k0 = 1;   % Near L4
Ytraj = Traj{k0};   % Nx4 [x y xdot ydot]
ttraj = Time{k0};

x0 = Ytraj(1,1);
y0 = Ytraj(1,2);
xd0 = Ytraj(1,3);
yd0 = Ytraj(1,4);

%% Display Lagrange info
disp('Lagrange points (x,y):');
for i = 1:5
    fprintf('L%d:  x = %12.8f,  y = %12.8f\n', i, Lx(i), Ly(i));
end

C0 = jacobi_const(x0,y0,xd0,yd0,mu);
fprintf('\nInitial Jacobi constant C0 = %.8f\n', C0);

%% Animation of chosen trajectory (rotating frame)
figure;
hold on;
plot(-mu,0,'ro','MarkerFaceColor','r');
plot(1-mu,0,'ko','MarkerFaceColor','k');

hTraj = plot(NaN,NaN,'b','LineWidth',1.2);
hP    = plot(NaN,NaN,'bo','MarkerFaceColor','b');

axis equal; grid on;
xlabel('x'); ylabel('y');
title(['Animation (rotating frame): ' labels{k0}]);

for n = 1:5:size(Ytraj,1)
    set(hTraj,'XData', Ytraj(1:n,1), 'YData', Ytraj(1:n,2));
    set(hP,   'XData', Ytraj(n,1),   'YData', Ytraj(n,2));
    drawnow;
end

%% Zero-velocity curve for chosen Jacobi constant
xs = linspace(-1.5,1.5,600);
ys = linspace(-1.2,1.2,400);
[Xg,Yg] = meshgrid(xs,ys);

r1 = sqrt((Xg+mu).^2 + Yg.^2);
r2 = sqrt((Xg-1+mu).^2 + Yg.^2);

Omega = 0.5*(Xg.^2 + Yg.^2) + (1-mu)./r1 + mu./r2;

figure;
contour(Xg, Yg, 2*Omega, [C0 C0], 'LineWidth',1.5);
hold on;
plot(-mu,0,'ro',1-mu,0,'ko');
axis equal; grid on;
xlabel('x'); ylabel('y');
title(['Zero-Velocity Curve, C = ' num2str(C0)]);

%% Poincaré section: y = 0 crossings
xcross = [];
xdotcross = [];

for n = 1:size(Ytraj,1)-1
    y1 = Ytraj(n,2);
    y2 = Ytraj(n+1,2);

    if y1*y2 < 0   % crosses y=0 (either direction)

        alpha = -y1 / (y2 - y1);

        x_at  = Ytraj(n,1) + alpha*(Ytraj(n+1,1) - Ytraj(n,1));
        xd_at = Ytraj(n,3) + alpha*(Ytraj(n+1,3) - Ytraj(n,3));

        xcross(end+1,1) = x_at;
        xdotcross(end+1,1) = xd_at;
    end
end

fprintf('Poincare points found: %d\n', length(xcross));

figure;
plot(xcross, xdotcross, 'k.', 'MarkerSize',8);
grid on;
xlabel('x');
ylabel('xdot');
title(['Poincare Section (y=0): ' labels{k0}]);

%% Functions
function [value,isterminal,direction] = yCrossEvent(~, y)
% yCrossEvent   ODE event function used to detect y = 0 crossings
% returns value = y. Not terminal; direction = 0 (both directions).
    value = y(2);        % y
    isterminal = 0;      % do not stop integration
    direction = 0;       % detect crossings in both directions
end

function [Lx, Ly] = lagrange_points(mu)
    fx = @(x) x ...
        - (1-mu)*(x+mu)./abs(x+mu).^3 ...
        - mu*(x-1+mu)./abs(x-1+mu).^3;

    L1x = fzero(fx, 0.7);
    L2x = fzero(fx, 1.2);
    L3x = fzero(fx, -1.0);

    xL4 = 0.5 - mu;
    yL4 =  sqrt(3)/2;
    xL5 = 0.5 - mu;
    yL5 = -sqrt(3)/2;

    Lx = [L1x; L2x; L3x; xL4; xL5];
    Ly = [0; 0; 0; yL4; yL5];
end

function dydt = cr3bp_ode(~, y, mu, perturb_handle, eps)
    % cr3bp_ode  Equations of motion for the planar CR3BP in rotating frame.
    %   dydt = cr3bp_ode(t, y, mu, perturb_handle, eps)
    %   y = [x; y; xdot; ydot]
    %
    %   perturb_handle is a function handle: [axp, ayp] = perturb_handle(x,y,xd,yd,mu,eps)
    %   if no perturbation is desired, set eps = 0 or pass [] for perturb_handle.

    x  = y(1); yy = y(2); u = y(3); v = y(4);
    % distances to primaries:
    r1 = sqrt((x + mu)^2 + yy^2);
    r2 = sqrt((x - 1 + mu)^2 + yy^2);

    % effective potential derivatives (Omega_x, Omega_y)
    Omega_x = x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
    Omega_y = yy - (1-mu)*yy/r1^3 - mu*yy/r2^3;

    % Equations of motion in rotating frame:
    xdd =  2*v + Omega_x; 
    ydd = -2*u + Omega_y;

    % apply optional small perturbation
    if nargin >= 5 && eps ~= 0 && ~isempty(perturb_handle)
        [axp, ayp] = perturb_handle(x, yy, u, v, mu, eps);
        xdd = xdd + axp;
        ydd = ydd + ayp;
    end

    dydt = [u; v; xdd; ydd];
end

function C = jacobi_const(x,y,xd,yd,mu)
% jacobi_const  Compute Jacobi constant (vectorized).
%   C = jacobi_const(x,y,xd,yd,mu)
%
%   x,y,xd,yd may be scalars or column vectors of equal length.


% Jacobi const is sort of like energy conservation. It is constant along
% certain trajectories. A higher C means less kinetic energy available. 

    r1 = sqrt((x + mu).^2 + y.^2); % The distance to the Earth
    r2 = sqrt((x - 1 + mu).^2 + y.^2); % Distance to the moon
    Omega = 0.5*(x.^2 + y.^2) + (1-mu)./r1 + mu./r2; % The effective potential
    C = 2.*Omega - (xd.^2 + yd.^2); % The Jacobi constant- denoting the zero-velocity curves where 2.*Omega = C
end

function [axp, ayp] = perturb_hook(x, y, xd, yd, mu, eps)
% perturb_hook  Example small perturbation hook.
%   [axp, ayp] = perturb_hook(x, y, xd, yd, mu, eps)
% Default example: small central 1/r^2 perturbation centered on the larger primary (-mu).
% Replace this function with your desired perturbation (radiation pressure, oblateness, thrust, etc).

    r1vec = [x + mu; y];
    r1 = norm(r1vec);
    if r1 == 0
        axp = 0; ayp = 0; return;
    end

    accel_mag = eps / r1^2;
    axp = -accel_mag * (r1vec(1)/r1);
    ayp = -accel_mag * (r1vec(2)/r1);
end

function [xcross, xdotcross] = poincare_section(Y)
% poincare_section  Find y=0 crossings with ydot > 0 (linear interpolation)
%   [xcross, xdotcross] = poincare_section(Y)
%   Y is Nx4 array: [x, y, xdot, ydot]
%
%   returns arrays of x and xdot at the crossing points.

    x = Y(:,1); y = Y(:,2); xd = Y(:,3); yd = Y(:,4);
    xcross = [];
    xdotcross = [];

    for k=1:length(y)-1
        if (y(k) <= 0 && y(k+1) > 0) || (y(k) < 0 && y(k+1) >= 0) || (y(k) >= 0 && y(k+1) < 0)
            % linear interpolation to find y=0 between k and k+1
            denom = (y(k+1) - y(k));
            if denom == 0
                continue;
            end
            tfrac = -y(k) / denom;
            x_at_cross = x(k) + tfrac*(x(k+1)-x(k));
            xd_at_cross = xd(k) + tfrac*(xd(k+1)-xd(k));
            yd_at_cross = yd(k) + tfrac*(yd(k+1)-yd(k));

            if yd_at_cross > 0
                xcross(end+1,1) = x_at_cross; %#ok<AGROW>
                xdotcross(end+1,1) = xd_at_cross; %#ok<AGROW>
            end
        end
    end
end