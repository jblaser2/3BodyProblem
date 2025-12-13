% animate_each_trajectory.m
% Animate each CR3BP trajectory in its own figure.

clear; close all; clc;

%% PARAMETERS
mu = 0.0121505856;      % Earth-Moon CR3BP
perturb_eps = 0.0;
tspan = [0, 200];
opts = odeset('RelTol',1e-9,'AbsTol',1e-11);

%% Lagrange points for convenient ICs
[Lx, Ly] = lagrange_points(mu);

%% Four initial conditions (rotating frame)
ic{1} = [Lx(1)-0.02; 0; 0; 0.03];        % near L1
ic{2} = [Lx(2)+0.02; 0; 0; -0.02];       % near L2
ic{3} = [Lx(4)+0.02; Ly(4)-0.02; 0; 0];  % near L4
ic{4} = [0.2; 0.05; -0.02; 0.02];        % transferring orbit

Ntraj = length(ic);

%% Time for synchronized animation
t_anim = linspace(tspan(1), tspan(2), 3500)';

%% Solve & animate each trajectory separately
for k = 1:Ntraj

    fprintf("Integrating trajectory %d...\n",k);
    [t, Y] = ode45(@(tt,yy) cr3bp_ode(tt,yy,mu,@perturb_hook,perturb_eps), ...
                    tspan, ic{k}, opts);

    % Interpolate onto common grid
    Yk = interp1(t, Y, t_anim, 'linear', 'extrap');

    xr = Yk(:,1); 
    yr = Yk(:,2);

    %% Convert rotating-frame -> inertial-frame coordinates
    theta = t_anim;
    X = xr .* cos(theta) - yr .* sin(theta);
    Yp = xr .* sin(theta) + yr .* cos(theta);

    %% Primaries orbiting in inertial frame
    prim1 = [-mu, 0];
    prim2 = [1-mu, 0];

    P1 = [prim1(1)*cos(theta) - prim1(2)*sin(theta), ...
          prim1(1)*sin(theta) + prim1(2)*cos(theta)];

    P2 = [prim2(1)*cos(theta) - prim2(2)*sin(theta), ...
          prim2(1)*sin(theta) + prim2(2)*cos(theta)];

    %% ---- NEW: One Figure Per Trajectory ----
    figure('Name',sprintf("Trajectory %d",k),'Color','w','NumberTitle','off',...
           'Position',[200 80 900 650]);
    ax = axes('NextPlot','add');
    axis equal; grid on;
    xlabel('X (inertial)'); ylabel('Y (inertial)');
    xlim([-1.5 1.5]); ylim([-1.2 1.2]);
    title(sprintf('CR3BP Trajectory %d (inertial frame)',k));

    % Trails & markers
    hTrail = plot(NaN,NaN,'b-','LineWidth',1.2);
    hObj   = plot(NaN,NaN,'bo','MarkerFaceColor','b','MarkerSize',6);
    hP1    = plot(NaN,NaN,'ro','MarkerFaceColor','r','MarkerSize',10);
    hP2    = plot(NaN,NaN,'ko','MarkerFaceColor','k','MarkerSize',8);

    trailLen = 300;

    fprintf("Animating trajectory %d...\n",k)
    for i = 1:length(t_anim)

        % primaries
        set(hP1,'XData',P1(i,1),'YData',P1(i,2));
        set(hP2,'XData',P2(i,1),'YData',P2(i,2));

        % trajectory
        idx0 = max(1, i-trailLen+1);
        set(hTrail,'XData',X(idx0:i),'YData',Yp(idx0:i));
        set(hObj,  'XData',X(i),     'YData',Yp(i));

        title(sprintf('Trajectory %d   t = %.2f',k,t_anim(i)));

        drawnow limitrate; 
    end
end

fprintf("Done. Each trajectory has its own figure.\n");
