function Kr = yaw_damper_gain(qbar)
%YAW_DAMPER_GAIN  Yaw damper gain scheduled on dynamic pressure.
%
%   Kr = yaw_damper_gain(qbar)      qbar = 0.5*rho*Va^2  [Pa]
%
%   The table is the output of design_gain_schedule() (target zeta 0.62,
%   washout a = 0.2 rad/s, rudder actuator tau_R included). Values are
%   linearly interpolated between breakpoints and held constant beyond
%   the ends (no extrapolation).
%
%   Shape: at low qbar (high alpha) the gain drops to ZERO. There the
%   directional stability is weak (Cn ~ 1 - alpha/15 deg), the Dutch roll
%   becomes roll dominated, and yaw rate feedback cannot add damping; any
%   gain only makes it worse. The gain peaks at moderate qbar and falls
%   again at high qbar as the rudder becomes more effective.

qbar_bp = [2210 2590 2990 3540 4100 4960 6020 8100 10760];
Kr_bp   = [0.00 0.58 0.91 1.08 1.12 1.06 0.90 0.76 0.67];

q  = min(max(qbar, qbar_bp(1)), qbar_bp(end));
Kr = interp1(qbar_bp, Kr_bp, q, 'linear');
end
