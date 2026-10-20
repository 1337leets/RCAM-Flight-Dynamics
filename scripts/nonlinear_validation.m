function R = nonlinear_validation(outdir)
%NONLINEAR_VALIDATION  Validate the yaw damper on the nonlinear 6-DOF model.
%
%   R = nonlinear_validation             % draws the figures, returns results
%   R = nonlinear_validation(outdir)     % also writes the figures to outdir
%
%   1) Rudder doublet (5 deg, 1 s + 1 s) at four flight points with the
%      damper off, with the fixed Kr = 1, and with the scheduled Kr(qbar).
%      The nonlinear response is compared with the response of the linear
%      lateral model to the same input.
%   2) The same test with a 10 deg doublet (large amplitude, actuator
%      saturation).
%   3) Turn: a 3 deg aileron step, to check whether yaw rate feedback with
%      and without washout opposes the turn.
%
%   The damper is the same as in yaw_damper_demo: washout a = 0.2 rad/s,
%   rudder actuator 1/(tau_R*s+1), rudder limit uR_lim. Elevator and
%   throttle are held at their trim values. Integration: fixed-step RK4.

if nargin < 1, outdir = ''; end
pts = [0 85; 0 120; 3000 70; 6000 140];          % [altitude m, Va m/s]
rho = @(h) interp1([0 3000 6000], [1.225 0.9093 0.6601], h);
dt = 0.01;  Tend = 30;
a  = 0.2;

R = struct([]);
for amp = [5 10]
    for k = 1:size(pts,1)
        p = rcam_params();  p.rho = rho(pts(k,1));  Va = pts(k,2);
        [~, x_tr, u_tr] = evalc('rcam_trim(Va, p)');
        [~,~,Alat,Blat] = rcam_linearize(x_tr, u_tr, p);
        din = @(t) deg2rad(amp) * ((t >= 1 & t < 2) - (t >= 2 & t < 3));

        modes = {'off', 'fixed', 'scheduled'};
        for m = 1:3
            ctl = struct('mode', modes{m}, 'a', a, 'washout', true);
            nl = sim_nonlinear(p, x_tr, u_tr, ctl, din, @(t) 0, dt, Tend);
            li = sim_linear(p, Alat, Blat, x_tr, ctl, din, dt, Tend);
            i  = numel(R) + 1;
            R(i).amp = amp;  R(i).alt = pts(k,1);  R(i).Va = Va;  R(i).mode = modes{m}; %#ok<AGROW>
            R(i).t = nl.t;  R(i).r_nl = nl.r;  R(i).r_li = li.r;
            R(i).dR = nl.dR;  R(i).beta = nl.beta;  R(i).Kr = nl.Kr;
            R(i).nrmse = norm(nl.r - li.r) / max(norm(li.r), eps);
            R(i).Jr    = trapz(nl.t, nl.r.^2);
            R(i).Jb    = trapz(nl.t, nl.beta.^2);
            R(i).sat = mean(abs(nl.dRcmd) >= p.uR_lim(2) - 1e-9);
        end
    end
end

% ---- Summary -------------------------------------------------------------
% Jr, Jb: int r^2 dt and int beta^2 dt, relative to the open loop at the
% same point. (Peak-to-peak logarithmic decrement is not used here: the r
% response is a sum of several modes and the damper adds a new damped
% pair, so the closed-loop eigenvalue is the right measure of Dutch roll
% damping.)
fprintf('\n%4s %6s %4s %-9s %6s | %8s | %7s %7s | %6s | %6s\n', 'amp','h','Va','mode','Kr', ...
        'NRMSE', 'Jr/Jr0', 'Jb/Jb0', 'dRmax', 'sat');
for i = 1:numel(R)
    i0 = find([R.amp] == R(i).amp & [R.alt] == R(i).alt & [R.Va] == R(i).Va & strcmp({R.mode}, 'off'));
    R(i).Jr_rel = R(i).Jr / R(i0).Jr;  R(i).Jb_rel = R(i).Jb / R(i0).Jb;
    fprintf('%4d %6d %4d %-9s %6.2f | %8.4f | %7.3f %7.3f | %6.1f | %5.1f%%\n', R(i).amp, R(i).alt, R(i).Va, ...
        R(i).mode, R(i).Kr(1), R(i).nrmse, R(i).Jr_rel, R(i).Jb_rel, max(abs(rad2deg(R(i).dR))), 100*R(i).sat);
end

T = turn_test(a);

f1 = plot_doublets(R, pts, 5);
f2 = plot_doublets(R, pts, 10);
f3 = plot_turn(T);
if ~isempty(outdir)
    exportgraphics(f1, fullfile(outdir, 'nl_doublet_5deg.png'),  'Resolution', 200);
    exportgraphics(f2, fullfile(outdir, 'nl_doublet_10deg.png'), 'Resolution', 200);
    exportgraphics(f3, fullfile(outdir, 'nl_turn.png'),          'Resolution', 200);
end
end

% =========================================================================
function S = sim_nonlinear(p, x_tr, u_tr, ctl, din, ain, dt, Tend)
% State z = [x (9) ; dR (actuator output, offset from trim) ; xw (washout)]
N = round(Tend/dt) + 1;  t = (0:N-1)'*dt;
z = [x_tr; 0; 0];
S.t = t;  S.r = zeros(N,1);  S.dR = zeros(N,1);  S.dRcmd = zeros(N,1);
S.beta = zeros(N,1);  S.phi = zeros(N,1);  S.Kr = zeros(N,1);
f = @(tt, zz) rhs(tt, zz, p, u_tr, ctl, din, ain);
for k = 1:N
    [~, out] = f(t(k), z);
    S.r(k) = z(6);  S.dR(k) = z(10);  S.dRcmd(k) = out.cmd;  S.Kr(k) = out.Kr;
    S.beta(k) = asin(z(2)/norm(z(1:3)));  S.phi(k) = z(7);
    if k == N, break; end
    k1 = f(t(k),        z);
    k2 = f(t(k) + dt/2, z + dt/2*k1);
    k3 = f(t(k) + dt/2, z + dt/2*k2);
    k4 = f(t(k) + dt,   z + dt*k3);
    z  = z + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end
end

function [zd, out] = rhs(t, z, p, u_tr, ctl, din, ain)
x = z(1:9);  dR = z(10);  xw = z(11);
r = x(6);                                   % r = 0 at trim
Va = norm(x(1:3));
Kr = gain(ctl.mode, 0.5*p.rho*Va^2);
if ctl.washout, y = r - ctl.a*xw; else, y = r; end
cmd = Kr*y + din(t);
cmd = max(min(cmd, p.uR_lim(2)), p.uR_lim(1));
u = u_tr;  u(1) = u_tr(1) + ain(t);  u(3) = u_tr(3) + dR;
zd = [rcam_dynamics(x, u, p);
      (cmd - dR)/p.tau_R;
      r - ctl.a*xw];
out.cmd = cmd;  out.Kr = Kr;
end

function S = sim_linear(p, Alat, Blat, x_tr, ctl, din, dt, Tend)
% Linear lateral model with the same damper; state [v p r phi | dR | xw]
N = round(Tend/dt) + 1;  t = (0:N-1)'*dt;
Kr = gain(ctl.mode, 0.5*p.rho*norm(x_tr(1:3))^2);   % fixed at the trim qbar
z = zeros(6,1);  S.r = zeros(N,1);
f = @(tt, zz) [Alat*zz(1:4) + Blat(:,2)*zz(5);
               (sat(Kr*(zz(3) - ctl.a*zz(6)) + din(tt), p) - zz(5))/p.tau_R;
               zz(3) - ctl.a*zz(6)];
for k = 1:N
    S.r(k) = z(3);
    if k == N, break; end
    k1 = f(t(k),        z);
    k2 = f(t(k) + dt/2, z + dt/2*k1);
    k3 = f(t(k) + dt/2, z + dt/2*k2);
    k4 = f(t(k) + dt,   z + dt*k3);
    z  = z + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end
end

function y = sat(u, p), y = max(min(u, p.uR_lim(2)), p.uR_lim(1)); end

function Kr = gain(mode, qbar)
switch mode
    case 'off',       Kr = 0;
    case 'fixed',     Kr = 1.0;                  % yaw_damper_demo
    case 'scheduled', Kr = yaw_damper_gain(qbar);
end
end

% =========================================================================
function T = turn_test(a)
% 0 m, 85 m/s: 3 deg aileron step (t >= 1 s), 40 s
p = rcam_params();
[~, x_tr, u_tr] = evalc('rcam_trim(85, p)');
ain = @(t) deg2rad(3) * (t >= 1);
cases = {'Damper OFF',                struct('mode','off',       'a',a, 'washout',true);
         'K_r(q) with washout',       struct('mode','scheduled', 'a',a, 'washout',true);
         'K_r(q) without washout',    struct('mode','scheduled', 'a',a, 'washout',false)};
T = struct([]);
for c = 1:size(cases,1)
    S = sim_nonlinear(p, x_tr, u_tr, cases{c,2}, @(t) 0, ain, 0.01, 40);
    T(c).name = cases{c,1};  T(c).t = S.t;  T(c).phi = S.phi;  T(c).r = S.r; %#ok<AGROW>
    T(c).dR = S.dR;  T(c).beta = S.beta;
end
fprintf('\nTurn test (0 m, 85 m/s, 3 deg aileron step), t = 40 s:\n');
for c = 1:numel(T)
    fprintf('  %-24s phi=%6.1f deg  r=%5.2f deg/s  dR=%5.2f deg  beta=%5.2f deg\n', T(c).name, ...
        rad2deg(T(c).phi(end)), rad2deg(T(c).r(end)), rad2deg(T(c).dR(end)), rad2deg(T(c).beta(end)));
end
end

% =========================================================================
function f = plot_doublets(R, pts, amp)
f = figure('Color','w','Position',[100 100 1150 720]);
if exist('theme','file'), theme(f,'light'); end    % R2025a+ defaults to dark
col = struct('off',[0.55 0.55 0.55], 'fixed',[0.85 0.33 0.10], 'scheduled',[0 0.45 0.74]);
for k = 1:size(pts,1)
    ax = subplot(2,2,k); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    for m = {'off','fixed','scheduled'}
        i = find([R.amp] == amp & [R.alt] == pts(k,1) & [R.Va] == pts(k,2) & strcmp({R.mode}, m{1}));
        plot(ax, R(i).t, rad2deg(R(i).r_nl), '-',  'Color', col.(m{1}), 'LineWidth', 2.0);
        plot(ax, R(i).t, rad2deg(R(i).r_li), '--', 'Color', col.(m{1}), 'LineWidth', 1.2, 'HandleVisibility','off');
    end
    ic = find([R.amp] == amp & [R.alt] == pts(k,1) & [R.Va] == pts(k,2) & strcmp({R.mode}, 'scheduled'));
    title(ax, sprintf('%d m, %d m/s  (K_r(q) = %.2f)', pts(k,1), pts(k,2), R(ic).Kr(1)));
    xlabel(ax, 'Time [s]');  ylabel(ax, 'r [\circ/s]');  set(ax, 'FontSize', 11);
    if k == 1, legend(ax, {'Damper OFF','Fixed K_r = 1','K_r(q)'}, 'Location','northeast'); end
end
sgtitle(sprintf(['%d\\circ rudder doublet  |  solid: nonlinear 6-DOF, ' ...
                 'dashed: linear lateral model'], amp), 'FontSize', 13, 'FontWeight', 'bold');
end

function f = plot_turn(T)
f = figure('Color','w','Position',[100 100 1150 650]);
if exist('theme','file'), theme(f,'light'); end
col = {[0.55 0.55 0.55], [0 0.45 0.74], [0.85 0.33 0.10]};
fld = {'phi','r','dR','beta'};
lab = {'Bank \phi [\circ]','Yaw rate r [\circ/s]','Rudder \delta_R [\circ]','Sideslip \beta [\circ]'};
for j = 1:4
    ax = subplot(2,2,j); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    for c = 1:numel(T)
        plot(ax, T(c).t, rad2deg(T(c).(fld{j})), '-', 'Color', col{c}, 'LineWidth', 2.0);
    end
    xlabel(ax, 'Time [s]');  ylabel(ax, lab{j});  set(ax, 'FontSize', 11);
    if j == 1, legend(ax, {T.name}, 'Location','northeast'); end
end
sgtitle('Turn test: 3\circ aileron step, 0 m, 85 m/s (nonlinear 6-DOF)', ...
        'FontSize', 13, 'FontWeight', 'bold');
end
