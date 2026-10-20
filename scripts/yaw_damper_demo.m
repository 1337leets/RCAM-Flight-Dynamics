function yaw_damper_demo(outfile)
%YAW_DAMPER_DEMO  Dutch roll with the yaw damper OFF and ON.
%
%   yaw_damper_demo              % writes the figure to yaw_damper_slide.png
%   yaw_damper_demo(outfile)     % writes the figure to the given file
%
%   Yaw damper law:  rudder_cmd = Kr * washout(yaw_rate),
%   with washout = s/(s+a). The rudder actuator is modeled as
%   1/(tau_R*s+1). Kr > 0 is the damping sign, because the rudder to yaw
%   rate entry of Blat is negative.
%
%   Kr and a were chosen from the closed-loop eigenvalues (actuator
%   included) by sweeping Kr = 0.25 to 4 and a = 0.05 to 1 at Va = 85 m/s.
%   With Kr = 1 and a = 0.2 the Dutch roll damping rises from 0.34 to
%   0.62. A smaller a gives more damping, but a slow washout also opposes
%   coordinated turns.

if nargin < 1, outfile = 'yaw_damper_slide.png'; end

p = rcam_params();
[x_tr, u_tr] = rcam_trim(85, p);
[~,~,Alat,Blat] = rcam_linearize(x_tr, u_tr, p);

IR = 3;                                   % yaw rate index (v,p,r,phi)

% ---- Isolate the Dutch roll (scaled to a yaw rate of 3 deg/s) ----------
[V,D] = eig(Alat); lam = diag(D);
cidx  = find(imag(lam) > 1e-6, 1);
vDR   = real(V(:,cidx));  vDR = vDR / vDR(IR);
x0    = vDR * deg2rad(3);

dt = 0.01; Tend = 40; N = round(Tend/dt); t = (0:N-1)*dt;

Kr = 1.0;                                 % positive means damping
a  = 0.2;                                 % washout corner [rad/s] (tau_w = 5 s)
tR = p.tau_R;                             % rudder actuator time constant [s]

% ---- Closed-loop eigenvalues: state [v p r phi | dR | xw] --------------
C   = [0 0 1 0];
Acl = [Alat,       Blat(:,2),  zeros(4,1);
       Kr*C/tR,    -1/tR,      -Kr*a/tR;
       C,          0,          -a];
z_off = dr_zeta(eig(Alat), lam(cidx));
z_on  = dr_zeta(eig(Acl),  lam(cidx));

% ---- 1) DAMPER OFF -----------------------------------------------------
x = x0; r_off = zeros(N,1);
for k = 1:N, r_off(k)=x(IR); x = x + dt*(Alat*x); end

% ---- 2) DAMPER ON ------------------------------------------------------
x = x0; xw = 0; dRa = 0; r_on = zeros(N,1); dR = zeros(N,1);
for k = 1:N
    r_on(k) = x(IR);
    yout = x(IR) - a*xw;                   % washout output
    xw   = xw + dt*(-a*xw + x(IR));
    uc   = Kr * yout;                      % rudder command
    uc   = max(min(uc, p.uR_lim(2)), p.uR_lim(1));
    dRa  = dRa + dt*(uc - dRa)/tR;         % actuator lag
    dR(k)= dRa;
    x    = x + dt*(Alat*x + Blat*[0; dRa]);
end

% ================= FIGURE ==============================================
f = figure('Color','w','Position',[100 100 1000 400]);
if exist('theme','file'), theme(f,'light'); end   % R2025a+ defaults to dark

ax1 = subplot(1,2,1);
plot(t, rad2deg(r_off), 'r-', 'LineWidth', 2.2); hold on;
plot(t, rad2deg(r_on),  'b-', 'LineWidth', 2.2);
grid on; box on;
xlabel('Time [s]','FontSize',12,'Color','k');
ylabel('Yaw rate  r  [\circ/s]','FontSize',12,'Color','k');
title(sprintf('Dutch Roll Damping  (\\zeta: %.2f \\rightarrow %.2f)', z_off, z_on), ...
      'FontSize',13,'FontWeight','bold','Color','k');
lg = legend({'Yaw damper OFF','Yaw damper ON'},'FontSize',11);
set(lg,'TextColor','k','Color','w');
set(ax1,'FontSize',11,'LineWidth',1,'XColor','k','YColor','k','GridColor',[.4 .4 .4]);

ax2 = subplot(1,2,2);
plot(t, rad2deg(dR), 'b-', 'LineWidth', 2);
grid on; box on;
xlabel('Time [s]','FontSize',12,'Color','k');
ylabel('Rudder deflection  \delta_R  [\circ]','FontSize',12,'Color','k');
title('Damper Rudder Activity','FontSize',13,'FontWeight','bold','Color','k');
set(ax2,'FontSize',11,'LineWidth',1,'XColor','k','YColor','k','GridColor',[.4 .4 .4]);

sgtitle('RCAM Lateral-Directional Model  |  V_a = 85 m/s', ...
        'FontSize',14,'FontWeight','bold','Color','k');

exportgraphics(f, outfile, 'Resolution', 200);

% ---- Summary -----------------------------------------------------------
fprintf('\n--- Yaw damper demo (Kr=%.2f, a=%.2f rad/s, tau_R=%.2f s) ---\n', Kr, a, tR);
fprintf('Dutch roll zeta: OFF = %.3f   ON = %.3f\n', z_off, z_on);
fprintf('Closed loop stable: %d\n', all(real(eig(Acl)) < 0));
fprintf('Peak rudder deflection: %.1f deg\n', max(abs(rad2deg(dR))));
fprintf('PNG: %s\n\n', outfile);
end

% -------------------------------------------------------------------------
function z = dr_zeta(lam, lam_dr0)
% Damping ratio of the complex pole closest to the open-loop Dutch roll
c = lam(imag(lam) > 1e-6);
[~, i] = min(abs(c - lam_dr0));
z = -real(c(i)) / abs(c(i));
end
