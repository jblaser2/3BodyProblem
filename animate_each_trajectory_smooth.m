% animate_each_trajectory_smooth.m
% Smooth, slowed-down per-trajectory CR3BP animation with fixed FPS and optional MP4 save.
% Requires: cr3bp_ode.m, lagrange_points.m, perturb_hook.m, jacobi_const.m

clear; close all; clc;

%% ---------- USER CONTROLS ----------
mu = 0.0121505856;
perturb_eps = 0.0;
tspan = [0, 200];        % nondim time
solver_opts = odeset('RelTol',1e-9,'AbsTol',1e-11);

fps = 24;                % frames per second for animation (try 20-30)
playback_speed = 1.0;    % 1.0 = real time (t units per second), <1 slower, >1 faster
save_video = false;      % true -> save an mp4 (can slow things down)
video_filename = 'traj_animation.mp4';

trail_seconds = 6;       % how many seconds of trail to show (in simulation time units)
max_trail_points = 400;  % cap to avoid huge lines (performance)

%% ---------- Initial conditions ----------
[Lx, Ly] = lagrange_points(mu);
ic{1} = [Lx(1)-0.02; 0; 0; 0.03];        % near L1
ic{2} = [Lx(2)+0.02; 0; 0; -0.02];       % near L2
ic{3} = [Lx(4)+0.02; Ly(4)-0.02; 0; 0];  % near L4
ic{4} = [0.2; 0.05; -0.02; 0.02];        % transfer-like
Ntraj = length(ic);

%% ---------- Integrate trajectories ----------
fprintf('Integrating %d trajectories...\n', Ntraj);
solT = cell(1,Ntraj);
solY = cell(1,Ntraj);
for k=1:Ntraj
    [t, Y] = ode45(@(tt,yy) cr3bp_ode(tt,yy,mu,@perturb_hook,perturb_eps), tspan, ic{k}, solver_opts);
    solT{k} = t(:);
    solY{k} = Y;
    fprintf('  traj %d: %d steps, t_end=%.3f\n', k, numel(t), t(end));
end

%% ---------- Create fixed frame timings ----------
sim_duration = tspan(2) - tspan(1);
frame_dt = 1 / fps;                      % seconds per frame (wall clock)
sim_time_per_frame = playback_speed * frame_dt; % simulation time units to advance per frame
Nframes = ceil(sim_duration / sim_time_per_frame);
t_frames = linspace(tspan(1), tspan(2), Nframes);

fprintf('Animating at %d fps, total frames %d, sim_time/frame %.4f\n', fps, Nframes, sim_time_per_frame);

%% ---------- Interpolate solutions onto frame times (use pchip for smoothness) ----------
Y_frames = NaN(Nframes, 4, Ntraj);
for k=1:Ntraj
    % Use PCHIP (shape-preserving cubic) for smooth motion and stable derivatives
    Y_frames(:,:,k) = interp1(solT{k}, solY{k}, t_frames, 'pchip');
end

%% ---------- Convert rotating frame -> inertial frame ----------
theta = t_frames;            % rotation angle = t (nondimensional omega=1)
cosT = cos(theta); sinT = sin(theta);

prim1_rot = [-mu, 0];
prim2_rot = [1-mu, 0];

P1_inert = [ prim1_rot(1)*cosT - prim1_rot(2)*sinT, ...
             prim1_rot(1)*sinT + prim1_rot(2)*cosT ];
P2_inert = [ prim2_rot(1)*cosT - prim2_rot(2)*sinT, ...
             prim2_rot(1)*sinT + prim2_rot(2)*cosT ];

Traj_inert = NaN(Nframes, 2, Ntraj);
for k=1:Ntraj
    xr = Y_frames(:,1,k); yr = Y_frames(:,2,k);
    Traj_inert(:,1,k) = xr .* cosT - yr .* sinT;   % X
    Traj_inert(:,2,k) = xr .* sinT + yr .* cosT;   % Y
end

%% ---------- Prepare video writer if requested ----------
if save_video
    vw = VideoWriter(video_filename,'MPEG-4');
    vw.FrameRate = fps;
    open(vw);
    fprintf('Recording to %s at %d fps...\n', video_filename, fps);
else
    vw = [];
end

%% ---------- Animation: one figure per trajectory (smooth update) ----------
for k = 1:Ntraj
    fig = figure('Name',sprintf('Trajectory %d (smooth)',k),'NumberTitle','off','Color','w', ...
                 'Position',[200 80 800 700]);
    ax = axes(fig,'NextPlot','add');
    axis equal; grid on;
    xlim([-1.5 1.5]); ylim([-1.2 1.2]);
    xlabel('X (inertial)'); ylabel('Y (inertial)');
    title(sprintf('CR3BP Trajectory %d (inertial frame)', k));

    % Preplot orbit guides and Lagrange markers
    thc = linspace(0,2*pi,400);
    plot(prim1_rot(1)*cos(thc), prim1_rot(1)*sin(thc), ':', 'LineWidth', 0.8);
    plot(prim2_rot(1)*cos(thc), prim2_rot(1)*sin(thc), ':', 'LineWidth', 0.8);
    Lx_all = [Lx(:)]; Ly_all = [Ly(:)];
    plot(Lx_all, Ly_all, 'g+','MarkerSize',8,'LineWidth',1.2);

    % animatedline for trail (limits to improve perf)
    trail_line = animatedline('LineWidth',1.5);
    % moving object marker
    hObj = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize', 8, 'MarkerEdgeColor','k');
    % primaries
    hP1 = plot(NaN, NaN, 'ro', 'MarkerFaceColor','r','MarkerSize',10);
    hP2 = plot(NaN, NaN, 'ko', 'MarkerFaceColor','k','MarkerSize',8);

    % compute trail length in frames (based on trail_seconds)
    trail_frames = min(max_trail_points, max(1, round(trail_seconds / sim_time_per_frame)));
    fprintf('Trajectory %d: trail frames = %d\n', k, trail_frames);

    % Animation loop for this one trajectory
    for i = 1:Nframes
        % current positions
        Xcur = Traj_inert(i,1,k);
        Ycur = Traj_inert(i,2,k);
        % primaries
        set(hP1, 'XData', P1_inert(i,1), 'YData', P1_inert(i,2));
        set(hP2, 'XData', P2_inert(i,1), 'YData', P2_inert(i,2));

        % update trail: using animatedline is faster & smoother than redrawing whole lines
        addpoints(trail_line, Xcur, Ycur);
        % trim points if too many
        % (animatedline has no direct trim, so we rebuild when length exceeds limit)
        if trail_line.NumPoints > trail_frames
            % extract current points, keep last trail_frames, recreate animatedline
            pts = trail_line.Points;          % Nx2 matrix
            pts = pts(end-trail_frames+1:end, :);
            clearpoints(trail_line);
            addpoints(trail_line, pts(:,1)', pts(:,2)');
        end

        % update marker
        set(hObj, 'XData', Xcur, 'YData', Ycur);

        % update title/time
        title(sprintf('Trajectory %d  t = %.2f', k, t_frames(i)));

        % capture or display frame
        drawnow limitrate;

        % write to video if enabled (grabframe)
        if save_video
            frame = getframe(fig);
            writeVideo(vw, frame);
        else
            % control real-time playback: pause for frame_dt seconds
            pause(frame_dt);   % can be removed if you want MATLAB to render as fast as possible
        end
    end

    % close figure if saving video and you want only one window
    % close(fig); 
end

if save_video
    close(vw);
    fprintf('Saved video: %s\n', video_filename);
end

fprintf('All animations finished.\n');
