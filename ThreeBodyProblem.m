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
eps = 0.01;

ICs = [ ...
    Lx(1)+eps,  Lx(2)-eps,  Lx(3)-eps,  Lx(4)+eps,  Lx(5)-eps;   % x
    Ly(1),      Ly(2),      Ly(3),      Ly(4)-eps,  Ly(5)+eps;   % y
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
xlabel('x'); ylabel('\dot{x}');
title('Phase space: x vs \dot{x}');
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
xlabel('y'); ylabel('\dot{y}');
title('Phase space: y vs \dot{y}');
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
    ylabel('\dot{x}');
    title(['Phase Space: x vs \dot{x}  (' labels{k} ')']);

    % mark initial condition
    hold on;
    plot(Traj{k}(1,1), Traj{k}(1,3), 'ro', 'MarkerFaceColor','r');
    legend('trajectory','initial condition','Location','best');

    % --- y vs ydot ---
    subplot(2,1,2);
    plot(Traj{k}(:,2), Traj{k}(:,4), 'LineWidth',1.4);
    grid on;
    xlabel('y');
    ylabel('\dot{y}');
    title(['Phase Space: y vs \dot{y}  (' labels{k} ')']);

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

k = 4;   % choose trajectory (e.g. near L4)

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

xlabel('$x$', 'Interpreter','latex');
ylabel('$\dot{x}$', 'Interpreter','latex');

title(['Poincaré Section ($y=0$, $\dot{y}>0$) -- ' labels{k}], ...
      'Interpreter','latex');

%% Choose ONE trajectory for detailed analysis
k0 = 4;   % Near L4
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
