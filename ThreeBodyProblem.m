% 3 Body Problem -- Modified

% Main script for the planar Circular Restricted Three-Body Problem (CR3BP)
% - calls modular functions: cr3bp_ode, lagrange_points, jacobi_const, ...
% - integrates, plots trajectory, Jacobi const, and Poincare section

clear; close all; clc;

%% ---------------- PARAMETERS ----------------
mu = 0.0121505856;    % mass ratio m2/(m1+m2) : Earth-Moon ~0.01215
perturb_eps = 0.0;    % perturbation amplitude (set >0 to enable)
tspan = [0, 200];     % nondimensional time (rotating frame units)
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'Events',@yCrossEvent);

%% ------------- initial condition --------------
[Lx, Ly] = lagrange_points(mu);
% Example: small offset from L4
x0 = Lx(4) + 0.01;
y0 = Ly(4) - 0.01;
xd0 = 0;
yd0 = 0;
y0_state = [x0; y0; xd0; yd0];

%% ------------- INTEGRATE ---------------------
[t, Y] = ode45(@(t,y) cr3bp_ode(t,y,mu,@perturb_hook,perturb_eps), tspan, y0_state, opts);

%% ------------- ANALYSIS ----------------------
C = jacobi_const(Y(:,1), Y(:,2), Y(:,3), Y(:,4), mu);

%% -------------- PLOTS -----------------
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

%% -------------- Poincaré section (y=0, ydot>0) ---------------
[px, pvx] = poincare_section(Y);
figure('Name','Poincare section (y=0, ydot>0)','NumberTitle','off');
plot(px, pvx, '.','MarkerSize',8);
xlabel('x'); ylabel('xdot'); title('Poincare section y=0 (ydot > 0)');
grid on;

%% ------------- Display Lagrange info --------------
disp('Lagrange points (x,y):');
for i=1:5
    fprintf('L%d:  x = %12.8f,  y = %12.8f\n', i, Lx(i), Ly(i));
end

fprintf('\nInitial Jacobi constant C0 = %.8f\n', jacobi_const(x0,y0,xd0,yd0,mu));
