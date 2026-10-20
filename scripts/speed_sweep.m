function T = speed_sweep(outfile)
%SPEED_SWEEP  Speed and altitude sweep of the lateral modes and the yaw damper.
%
%   T = speed_sweep              % returns the table and draws the figure
%   T = speed_sweep(outfile)     % also writes the figure to a file
%
%   At every (Va, rho) point the aircraft is trimmed for straight and
%   level flight and linearized. The yaw damper (washout a = 0.2, rudder
%   actuator tau_R) is closed in two ways:
%     - with the fixed gain Kr = 1 (yaw_damper_demo, tuned at 85 m/s)
%     - with the scheduled gain Kr(qbar) (yaw_damper_gain)
%
%   The RCAM atmosphere is constant, so altitude is represented by the ISA
%   density only. Thrust does not lapse with altitude in the RCAM engine
%   model.

if nargin < 1, outfile = ''; end

Va_list  = 60:5:140;                       % m/s
alt_list = [0 3000 6000];                  % m
rho_list = [1.225 0.9093 0.6601];          % kg/m^3, ISA
Kr = 1.0;  a = 0.2;                        % yaw damper (yaw_damper_demo)

p0 = rcam_params();
n  = numel(Va_list) * numel(alt_list);
T  = table('Size',[n 17], 'VariableTypes',[repmat({'double'},1,16), {'logical'}], ...
    'VariableNames',{'alt_m','Va','qbar','alpha_deg','uE_deg','uTH_deg', ...
                     'dr_wn','dr_zeta','dr_zeta_yd','dr_zeta_gs', ...
                     'spiral_lam','spiral_lam_yd','spiral_lam_gs', ...
                     'Kr_gs','roll_tau','residual','trim_ok'});
k = 0;
for ia = 1:numel(alt_list)
    p = p0;  p.rho = rho_list(ia);
    for Va = Va_list
        k = k + 1;
        [~, x_tr, u_tr, info] = evalc('rcam_trim(Va, p)');   % silence the trim printout
        [~,~,Alat,Blat] = rcam_linearize(x_tr, u_tr, p);

        lam = eig(Alat);
        dr  = lam(imag(lam) > 1e-6);
        re  = sort(real(lam(abs(imag(lam)) <= 1e-6)));   % [roll; spiral]
        qbar = 0.5*p.rho*Va^2;
        Kgs  = yaw_damper_gain(qbar);
        lcl  = eig(closed_loop(Alat, Blat(:,2), p.tau_R, Kr,  a));
        lgs  = eig(closed_loop(Alat, Blat(:,2), p.tau_R, Kgs, a));

        ok = info.residual < 1e-6 && ...
             u_tr(2) >= p.uE_lim(1)  && u_tr(2) <= p.uE_lim(2) && ...
             u_tr(4) >= p.uTH_lim(1) && u_tr(4) <= p.uTH_lim(2) && ...
             atan2(x_tr(3), x_tr(1)) <= p.alpha_switch && ...
             numel(dr) == 1 && numel(re) == 2;

        T(k,:) = {alt_list(ia), Va, qbar, info.alpha_deg, info.uE_deg, info.uTH_deg, ...
                  nan, nan, nan, nan, nan, nan, nan, Kgs, nan, info.residual, ok};
        if numel(dr) == 1 && numel(re) == 2
            T.dr_wn(k)         = abs(dr);
            T.dr_zeta(k)       = -real(dr)/abs(dr);
            T.dr_zeta_yd(k)    = dr_zeta(lcl, dr);
            T.dr_zeta_gs(k)    = dr_zeta(lgs, dr);
            T.spiral_lam(k)    = re(2);
            T.spiral_lam_yd(k) = slowest_real(lcl);
            T.spiral_lam_gs(k) = slowest_real(lgs);
            T.roll_tau(k)      = -1/re(1);
        end
    end
end

% ---- Summary -------------------------------------------------------------
fprintf('\n%6s %5s %6s %7s %7s | %6s %6s %6s %6s %5s | %8s %8s | %s\n', 'h[m]','Va','qbar','alpha','uE', ...
        'wn_DR','z_DR','z_YD','z_GS','Kr_GS','spiral','spir_GS','trim');
for k = 1:height(T)
    fprintf('%6d %5d %6.0f %7.2f %7.2f | %6.3f %6.3f %6.3f %6.3f %5.2f | %+8.4f %+8.4f | %s\n', ...
        T.alt_m(k), T.Va(k), T.qbar(k), T.alpha_deg(k), T.uE_deg(k), ...
        T.dr_wn(k), T.dr_zeta(k), T.dr_zeta_yd(k), T.dr_zeta_gs(k), T.Kr_gs(k), ...
        T.spiral_lam(k), T.spiral_lam_gs(k), string(ifelse(T.trim_ok(k), 'ok', 'OUT OF LIMITS')));
end

plot_sweep(T, alt_list, outfile);
end

% -------------------------------------------------------------------------
function s = ifelse(c, a, b), if c, s = a; else, s = b; end, end

function Acl = closed_loop(Alat, bR, tR, Kr, a)
% state [v p r phi | dR | xw],  dR_cmd = Kr*(r - a*xw)  (same as yaw_damper_demo)
C   = [0 0 1 0];
Acl = [Alat,     bR,      zeros(4,1);
       Kr*C/tR, -1/tR,   -Kr*a/tR;
       C,        0,      -a];
end

function z = dr_zeta(lam, lam_dr0)
% Damping ratio of the complex pole closest to the open-loop Dutch roll
c = lam(imag(lam) > 1e-6);
if isempty(c), z = 1; return; end
[~, i] = min(abs(c - lam_dr0));
z = -real(c(i)) / abs(c(i));
end

function s = slowest_real(lam)
% Real pole closest to the origin (the spiral mode in closed loop)
r = real(lam(abs(imag(lam)) <= 1e-6));
[~, i] = min(abs(r));
s = r(i);
end

function plot_sweep(T, alt_list, outfile)
f = figure('Color','w','Position',[100 100 1100 750]);
if exist('theme','file'), theme(f,'light'); end    % R2025a+ defaults to dark
c  = [0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];
lg = arrayfun(@(h) sprintf('%d m', h), alt_list, 'UniformOutput', false);

ax = gobjects(4,1);
for i = 1:4, ax(i) = subplot(2,2,i); hold(ax(i),'on'); grid(ax(i),'on'); box(ax(i),'on'); end
for ia = 1:numel(alt_list)
    S  = T(T.alt_m == alt_list(ia), :);
    vars = {'dr_zeta','dr_zeta_yd','dr_zeta_gs','dr_wn','spiral_lam','spiral_lam_yd', ...
            'spiral_lam_gs','Kr_gs','uE_deg','alpha_deg'};
    S{~S.trim_ok, vars} = NaN;              % do not plot points outside the trim limits
    plot(ax(1), S.Va, S.dr_zeta,    '--', 'Color', c(ia,:), 'LineWidth', 1.5, 'HandleVisibility','off');
    plot(ax(1), S.Va, S.dr_zeta_yd, ':',  'Color', c(ia,:), 'LineWidth', 1.8, 'HandleVisibility','off');
    plot(ax(1), S.Va, S.dr_zeta_gs, '-',  'Color', c(ia,:), 'LineWidth', 2.2);
    plot(ax(2), S.Va, S.Kr_gs,      '-',  'Color', c(ia,:), 'LineWidth', 2.2);
    plot(ax(3), S.Va, S.spiral_lam,    '--', 'Color', c(ia,:), 'LineWidth', 1.5, 'HandleVisibility','off');
    plot(ax(3), S.Va, S.spiral_lam_gs, '-',  'Color', c(ia,:), 'LineWidth', 2.2);
    plot(ax(4), S.Va, S.uE_deg,     '-',  'Color', c(ia,:), 'LineWidth', 2.2);
    plot(ax(4), S.Va, S.alpha_deg,  ':',  'Color', c(ia,:), 'LineWidth', 1.8, 'HandleVisibility','off');
end
yline(ax(1), 0.19, 'k-.', 'MIL-F-8785C L1 Cat. A: 0.19', 'LabelHorizontalAlignment','left', 'HandleVisibility','off');
yline(ax(1), 0.6,  'k:',  'target 0.6', 'LabelHorizontalAlignment','left', 'HandleVisibility','off');
yline(ax(2), 1.0,  'k:',  'fixed K_r = 1', 'LabelHorizontalAlignment','left', 'HandleVisibility','off');
yline(ax(3), 0,    'k-', 'HandleVisibility','off');
yline(ax(4), -25,  'k:', 'uE lower limit', 'LabelHorizontalAlignment','left', 'HandleVisibility','off');

title(ax(1), 'Dutch roll \zeta  (dashed: OFF, dotted: fixed K_r, solid: K_r(q))');
ylabel(ax(1), '\zeta_{DR}');
title(ax(2), 'Scheduled gain K_r(\itq\rm)');  ylabel(ax(2), 'K_r [rad/(rad/s)]');
title(ax(3), 'Spiral pole  (>0 is unstable; dashed: OFF, solid: K_r(q))');
ylabel(ax(3), '\lambda_{spiral} [1/s]');
title(ax(4), 'Trim: uE (solid), \alpha (dotted)');  ylabel(ax(4), '[deg]');
for i = 1:4, xlabel(ax(i), 'V_a [m/s]'); set(ax(i), 'FontSize', 11); end
legend(ax(2), lg, 'Location', 'southeast');
sgtitle(sprintf('RCAM speed/altitude sweep  |  yaw damper a = 0.2 rad/s, \\tau_R = 0.3 s  |  points outside the trim limits are not plotted'), ...
        'FontSize', 13, 'FontWeight', 'bold');
if ~isempty(outfile), exportgraphics(f, outfile, 'Resolution', 200); end
end
