function [x_tr, u_tr, info] = rcam_trim(Va, p)
%RCAM_TRIM  Find the trim point for straight and level flight.
%
%   [x_tr, u_tr, info] = rcam_trim(Va, p)
%
%   Trimming means finding the (x, u) pair for which xdot is close to
%   zero: forces and moments balance and no state changes. All of the
%   linear control analysis in this repository is built around that point.
%
%   With wings level and zero flight path angle, only three variables are
%   free:
%     alpha  so that lift balances weight
%     uE     so that the pitching moment is zero
%     uTH    so that thrust balances drag
%   Everything else follows from these (theta = alpha because gamma = 0).

if nargin < 2, p = rcam_params(); end

% Initial guess: [alpha; uE; uTH]  (rad)
z0 = [0.08; -0.10; 0.08];

cost = @(z) trim_cost(z, Va, p);
opts = optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',5e4,'MaxIter',5e4,...
                'Display','off');
[z, fval, flag] = fminsearch(cost, z0, opts);

[x_tr, u_tr] = build_state(z, Va);

info.residual = fval;
info.exitflag = flag;
info.alpha_deg = z(1)*180/pi;
info.uE_deg    = z(2)*180/pi;
info.uTH_deg   = z(3)*180/pi;
info.Va        = Va;

fprintf('--- TRIM @ Va = %.1f m/s ---\n', Va);
fprintf('  alpha = %7.3f deg\n', info.alpha_deg);
fprintf('  uE    = %7.3f deg\n', info.uE_deg);
fprintf('  uTH   = %7.3f deg (both engines)\n', info.uTH_deg);
fprintf('  residual = %.3e   %s\n', fval, ...
        ternary(fval < 1e-6, '[OK]', '[!! not converged]'));
end

% -------------------------------------------------------------------------
function J = trim_cost(z, Va, p)
[x, u] = build_state(z, Va);
xdot = rcam_dynamics(x, u, p);
% Only the first six states matter: velocities and rates must be constant.
% (The phi/theta/psi derivatives are zero by construction.)
J = sum(xdot(1:6).^2);
end

% -------------------------------------------------------------------------
function [x, u] = build_state(z, Va)
alpha = z(1);  uE = z(2);  uTH = z(3);

% gamma = 0 (level flight) => theta = alpha
theta = alpha;

x = [ Va*cos(alpha);   % u
      0;               % v   (no sideslip)
      Va*sin(alpha);   % w
      0; 0; 0;         % p q r
      0;               % phi (wings level)
      theta;
      0 ];             % psi
u = [0; uE; 0; uTH; uTH];   % zero aileron and rudder, symmetric throttle
end

% -------------------------------------------------------------------------
function out = ternary(c,a,b), if c, out=a; else, out=b; end, end
