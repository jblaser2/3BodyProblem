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
tspan = [0, 200];     % nondimensional time (rotating frame units)
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'Events',@yCrossEvent);

%% initial condition 
[Lx, Ly] = lagrange_points(mu);
% Example: small offset from L4
x0 = Lx(4) + 0.01;
y0 = Ly(4) - 0.01;
xd0 = 0;
yd0 = 0;
y0_state = [x0; y0; xd0; yd0];

%% Integrate
[t, Y] = ode45(@(t,y) cr3bp_ode(t,y,mu,@perturb_hook,perturb_eps), tspan, y0_state, opts);

%% Analysis
C = jacobi_const(Y(:,1), Y(:,2), Y(:,3), Y(:,4), mu);

%% Plots
figure('Name','CR3BP trajectory in rotating frame','NumberTitle','off','Position',[100 100 900 400]);
subplot(1,2,1);
plot(Y(:,1), Y(:,2), 'b'); hold on;
plot(-mu, 0, 'ro','MarkerFaceColor','r','DisplayName','Primary 1'); % m1
plot(1-mu, 0, 'ko','MarkerFaceColor','k','DisplayName','Primary 2'); % m2
plot(Lx, Ly, 'g+','MarkerSize',10,'LineWidth',1.5);
axis equal;
xlabel('x'); ylabel('y'); title('Trajectory (rotating frame)');
legend('trajectory','Primary1','Primary2','Lagrange points','Location','best');

subplot(1,2,2);
plot(t, C, '-'); xlabel('t'); ylabel('Jacobi constant C'); title('Jacobi constant along trajectory');

%%  Poincaré section (y=0, ydot>0)
[px, pvx] = poincare_section(Y);
figure('Name','Poincare section (y=0, ydot>0)','NumberTitle','off');
plot(px, pvx, '.','MarkerSize',8);
xlabel('x'); ylabel('xdot'); title('Poincare section y=0 (ydot > 0)');
grid on;

%% Display Lagrange info 
disp('Lagrange points (x,y):');
for i=1:5
    fprintf('L%d:  x = %12.8f,  y = %12.8f\n', i, Lx(i), Ly(i));
end

fprintf('\nInitial Jacobi constant C0 = %.8f\n', jacobi_const(x0,y0,xd0,yd0,mu));


%% Y is Nx4 from ode45: [x,y,xd,yd]
figure;
plot(-mu,0,'ro','MarkerFaceColor','r'); hold on;    % primary 1
plot(1-mu,0,'ko','MarkerFaceColor','k');            % primary 2
hTraj = plot(NaN,NaN,'b');                         % trajectory as it's drawn
hP = plot(NaN,NaN,'b.','MarkerSize',8);            % moving point
axis equal; xlabel('x'); ylabel('y'); title('Animation (rotating frame)');

for k=1:5:length(Y)   % step through frames (skip frames if dense)
    set(hTraj,'XData', Y(1:k,1), 'YData', Y(1:k,2));
    set(hP,'XData', Y(k,1), 'YData', Y(k,2));
    drawnow;
    pause(0.01);
end


%% choose grid
xs = linspace(-1.5,1.5,600);
ys = linspace(-1.0,1.0,400);
[X,Y] = meshgrid(xs,ys);
r1 = sqrt((X+mu).^2 + Y.^2);
r2 = sqrt((X-1+mu).^2 + Y.^2);
Omega = 0.5*(X.^2 + Y.^2) + (1-mu)./r1 + mu./r2;

C0 = jacobi_const(x0,y0,xd0,yd0,mu);   % pick initial state's Jacobi constant
% plot contour 2*Omega = C0
figure; contour(X,Y,2*Omega, [C0 C0], 'LineWidth',1.5);
hold on;
plot(-mu,0,'ro',1-mu,0,'ko'); axis equal;
xlabel('x'); ylabel('y'); title(['Zero-velocity curve, C = ' num2str(C0)]);

%% Y: Nx4, columns x,y,xd,yd
xcross = []; xdotcross = [];
for k=1:size(Y,1)-1
    yk = Y(k,2); yk1 = Y(k+1,2);
    if yk<=0 && yk1>0  % upward crossing
        tfrac = -yk/(yk1-yk);
        x_at = Y(k,1) + tfrac*(Y(k+1,1)-Y(k,1));
        xd_at = Y(k,3) + tfrac*(Y(k+1,3)-Y(k,3));
        xcross(end+1,1)=x_at;
        xdotcross(end+1,1)=xd_at;
    end
end
figure; plot(xcross, xdotcross, '.'); xlabel('x'); ylabel('xdot'); title('Poincare: y=0, ydot>0');

