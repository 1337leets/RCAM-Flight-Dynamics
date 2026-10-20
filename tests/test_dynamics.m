function test_dynamics()
%TEST_DYNAMICS  Physical consistency checks for rcam_dynamics.m.
%
%   These tests do not compare against a reference solution. They check
%   that the physics is self-consistent, and a failing test points to
%   where an error is likely to be.

p = rcam_params();
fprintf('\n=========== rcam_dynamics tests ===========\n');
np = 0; nf = 0;

% ---------------------------------------------------------------- TEST 1
% Output size and finiteness
x = [85; 0; 6; 0; 0; 0; 0; 0.07; 0];
u = [0; -0.1; 0; 0.08; 0.08];
try
    xdot = rcam_dynamics(x, u, p);
    ok = isequal(size(xdot), [9 1]) && all(isfinite(xdot));
catch ME
    ok = false;  fprintf('  error: %s\n', ME.message);
end
[np,nf] = report('1. Output is 9x1 and finite', ok, np, nf, ...
    'rcam_dynamics does not run, or it returns NaN.');

if nf > 0, summary(np,nf); return; end

% ---------------------------------------------------------------- TEST 2
% SYMMETRY: with beta = 0, phi = 0, zero aileron and rudder, and
% symmetric thrust, the lateral motion must not be excited, so vdot,
% pdot and rdot must be about zero.
% This test catches most sign and index errors.
x = [85; 0; 6; 0; 0; 0; 0; 0.07; 0];
u = [0; -0.1; 0; 0.08; 0.08];
xdot = rcam_dynamics(x, u, p);
lat = xdot([2 4 6 7]);          % vdot, pdot, rdot, phidot
ok = max(abs(lat)) < 1e-9;
[np,nf] = report('2. Symmetric flight does not excite the lateral axis', ok, np, nf, ...
    sprintf(['Lateral derivatives are not zero (max=%.3e). Likely causes:\n' ...
             '     - beta terms in the wrong place in the eta vector\n' ...
             '     - cross() arguments in the wrong order\n' ...
             '     - wrong signs in the engine moment arms (mew1/mew2)'], max(abs(lat))));

% ---------------------------------------------------------------- TEST 3
% ASYMMETRIC THRUST: with the left engine off, a yawing moment must
% appear. Its sign may be either way, but it must not be zero (it is
% zero when the engine moment arm is missing).
u_ae = [0; -0.1; 0; 0.0; 0.16];   % left engine off, right engine doubled
xdot_ae = rcam_dynamics(x, u_ae, p);
ok = abs(xdot_ae(6)) > 1e-4;      % rdot
[np,nf] = report('3. Asymmetric thrust produces a yawing moment', ok, np, nf, ...
    'rdot is about zero. The engine moment arm may be missing.');

% ---------------------------------------------------------------- TEST 4
% Does the trim converge?
Va = 85;
[x_tr, u_tr, info] = rcam_trim(Va, p);
ok = info.residual < 1e-6;
[np,nf] = report('4. Trim converges (residual < 1e-6)', ok, np, nf, ...
    'No trim point found, which points to a structural error in the dynamics.');

% ---------------------------------------------------------------- TEST 5
% Is the trim physically reasonable?
ok = info.alpha_deg > 0 && info.alpha_deg < 12 && ...
     info.uTH_deg  > 0 && info.uTH_deg  < 15;
[np,nf] = report('5. Trim values are physical', ok, np, nf, ...
    sprintf(['alpha=%.2f deg, uTH=%.2f deg is not plausible.\n' ...
             '     An airliner at this speed flies at an alpha of roughly 0 to 8 deg.\n' ...
             '     A negative alpha, or one above 20 deg, means the lift and weight\n' ...
             '     balance has the wrong sign.'], info.alpha_deg, info.uTH_deg));

% ---------------------------------------------------------------- TEST 6
% The lateral-directional A matrix must have exactly one complex pair
% (the Dutch roll)
[~,~,Alat] = rcam_linearize(x_tr, u_tr, p);
lam = eig(Alat);
ncomplex = sum(abs(imag(lam)) > 1e-6);
ok = (ncomplex == 2);
[np,nf] = report('6. Exactly one complex lateral pair (Dutch roll)', ok, np, nf, ...
    sprintf('Found %d complex eigenvalues, expected 2.', ncomplex));

% ---------------------------------------------------------------- TEST 7
% Is the Dutch roll frequency in a reasonable range?
if ncomplex == 2
    dr = lam(abs(imag(lam)) > 1e-6);
    wn = abs(dr(1));
    ok = wn > 0.3 && wn < 3.0;
    [np,nf] = report('7. Dutch roll frequency is plausible (0.3 to 3 rad/s)', ok, np, nf, ...
        sprintf('wn = %.3f rad/s is unusual for a transport aircraft.', wn));

    zeta = -real(dr(1))/abs(dr(1));
    fprintf('\n  >>> Dutch roll:  wn = %.3f rad/s,  zeta = %.4f\n', wn, zeta);
end

summary(np, nf);
end

% -------------------------------------------------------------------------
function [np,nf] = report(name, ok, np, nf, hint)
if ok
    fprintf('  [PASS]  %s\n', name);  np = np + 1;
else
    fprintf('  [FAIL]  %s\n', name);
    fprintf('     -> %s\n', hint);     nf = nf + 1;
end
end

function summary(np, nf)
fprintf('---------------------------------------------\n');
fprintf('  %d passed, %d failed\n', np, nf);
if nf == 0
    fprintf('  The model is consistent.\n');
end
fprintf('=============================================\n\n');
end
