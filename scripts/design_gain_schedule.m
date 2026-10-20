function [qbar_bp, Kr_bp, D] = design_gain_schedule(zeta_target)
%DESIGN_GAIN_SCHEDULE  Design the yaw damper gain schedule Kr(qbar).
%
%   [qbar_bp, Kr_bp, D] = design_gain_schedule          % zeta_target = 0.62
%   [qbar_bp, Kr_bp, D] = design_gain_schedule(zeta_target)
%
%   At every valid trim point (the same envelope as speed_sweep), Kr is
%   swept with the washout a = 0.2 and the rudder actuator tau_R fixed:
%     - if zeta_target can be reached, the SMALLEST Kr that reaches it is
%       chosen, so the rudder is not used more than necessary;
%     - otherwise the SMALLEST Kr that comes within 0.005 of the best
%       achievable zeta is chosen. At low qbar and high alpha the curve is
%       flat, so this rule fades the damper out on its own ("do no harm").
%   The default target of 0.62 is the desired 0.6 plus a margin for the
%   interpolation loss between breakpoints.
%   The points are binned by dynamic pressure (qbar), and the median of
%   each bin becomes a breakpoint. In straight and level flight
%   qbar = W/(CL*S), so qbar also fixes the angle of attack, which makes
%   it a good scheduling variable.
%
%   The resulting table is hard-coded in yaw_damper_gain.m.

if nargin < 1, zeta_target = 0.62; end

Va_list  = 60:5:140;
rho_list = [1.225 0.9093 0.6601];          % 0 / 3000 / 6000 m, ISA
a  = 0.2;
Ks = 0:0.02:3;

p0 = rcam_params();
D  = zeros(0,5);                            % [qbar Va rho Kr* zeta*]
for rho = rho_list
    p = p0;  p.rho = rho;
    for Va = Va_list
        [~, x_tr, u_tr, info] = evalc('rcam_trim(Va, p)');
        if ~trim_ok(x_tr, u_tr, info, p), continue; end
        [~,~,Alat,Blat] = rcam_linearize(x_tr, u_tr, p);
        lam = eig(Alat);  dr0 = lam(imag(lam) > 1e-6);
        z = arrayfun(@(K) yd_zeta(Alat, Blat(:,2), p.tau_R, K, a, dr0), Ks);
        j = find(z >= zeta_target, 1);
        if isempty(j), j = find(z >= max(z) - 0.005, 1); end
        D(end+1,:) = [0.5*rho*Va^2, Va, rho, Ks(j), z(j)]; %#ok<AGROW>
    end
end

% ---- qbar bins -> breakpoints --------------------------------------------
edges   = [0 2400 2800 3300 3800 4500 5500 7000 9000 inf];
qbar_bp = zeros(1, numel(edges)-1);  Kr_bp = qbar_bp;
for i = 1:numel(edges)-1
    in = D(:,1) >= edges(i) & D(:,1) < edges(i+1);
    qbar_bp(i) = round(median(D(in,1)), -1);
    Kr_bp(i)   = round(median(D(in,4)), 2);
end

fprintf('\nKr(qbar) table (zeta_target = %.2f, a = %.2f rad/s):\n', zeta_target, a);
fprintf('  qbar [Pa]: %s\n', sprintf('%7.0f', qbar_bp));
fprintf('  Kr       : %s\n', sprintf('%7.2f', Kr_bp));
end

% -------------------------------------------------------------------------
function ok = trim_ok(x_tr, u_tr, info, p)
ok = info.residual < 1e-6 && ...
     u_tr(2) >= p.uE_lim(1)  && u_tr(2) <= p.uE_lim(2) && ...
     u_tr(4) >= p.uTH_lim(1) && u_tr(4) <= p.uTH_lim(2) && ...
     atan2(x_tr(3), x_tr(1)) <= p.alpha_switch;
end

function z = yd_zeta(Alat, bR, tR, Kr, a, dr0)
% Closed loop [v p r phi | dR | xw], dR_cmd = Kr*(r - a*xw); -1 if unstable
C = [0 0 1 0];
A = [Alat, bR, zeros(4,1); Kr*C/tR, -1/tR, -Kr*a/tR; C, 0, -a];
l = eig(A);
if any(real(l) >= 0), z = -1; return; end
c = l(imag(l) > 1e-6);
if isempty(c), z = 1; return; end
[~, i] = min(abs(c - dr0));
z = -real(c(i)) / abs(c(i));
end
