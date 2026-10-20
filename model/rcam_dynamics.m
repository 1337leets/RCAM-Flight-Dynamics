function xdot = rcam_dynamics(x, u, p)
%RCAM_DYNAMICS  Nonlinear 6-DOF RCAM aircraft dynamics.
%
%   xdot = rcam_dynamics(x, u, p)
%
%   STATE VECTOR x (9x1), all in BODY axes:
%     x(1) = u      [m/s]    forward velocity
%     x(2) = v      [m/s]    lateral velocity
%     x(3) = w      [m/s]    vertical velocity
%     x(4) = p      [rad/s]  roll rate
%     x(5) = q      [rad/s]  pitch rate
%     x(6) = r      [rad/s]  yaw rate
%     x(7) = phi    [rad]    bank angle
%     x(8) = theta  [rad]    pitch angle
%     x(9) = psi    [rad]    heading
%
%   CONTROL VECTOR u (5x1):
%     u(1) = uA     [rad]    aileron
%     u(2) = uE     [rad]    elevator (an all-moving tailplane in RCAM)
%     u(3) = uR     [rad]    rudder
%     u(4) = uTH1   [rad]    left throttle
%     u(5) = uTH2   [rad]    right throttle
%
%   Limitations: the Euler angle kinematics are singular at theta = +/-90
%   deg, and the air data are undefined at Va = 0 (the result is NaN).
%
% ---- 0. Unpack ----------------------------------------------------------
u_b = x(1);  v_b = x(2);  w_b = x(3);
p_b = x(4);  q_b = x(5);  r_b = x(6);
phi = x(7);  theta = x(8);  % psi = x(9) does not enter the dynamics

uA = u(1);  uE = u(2);  uR = u(3);  uTH1 = u(4);  uTH2 = u(5);

V_b   = [u_b; v_b; w_b];       % velocity vector in body axes
wbe_b = [p_b; q_b; r_b];       % angular velocity vector in body axes


%% ---- 1. Air data -------------------------------------------------------
Va = norm(V_b);				   % airspeed
alpha = atan2(w_b, u_b);       % angle of attack
beta = asin(v_b / Va);         % sideslip angle
Q = 0.5 * p.rho * Va^2;        % dynamic pressure


%% ---- 2. Aerodynamic coefficients ---------------------------------------
if alpha <= p.alpha_switch
    CL_wb = p.n * (alpha - p.alpha_L0);
else
    CL_wb = p.a3*alpha^3 + p.a2*alpha^2 + p.a1*alpha + p.a0;
end

epsilon = p.depsda * (alpha - p.alpha_L0);
alpha_t = alpha - epsilon + uE + 1.3*q_b*p.lt/Va;
CL_t    = p.nt * (p.St/p.S) * alpha_t;

CL = CL_wb + CL_t;

CD = p.CDmin + p.d1 * (p.n*alpha + p.d0)^2;

CY = p.CY_beta * beta + p.CY_dR * uR;


%% ---- 3. Aerodynamic forces ---------------------------------------------
FA_s = [ -CD * Q * p.S ;
		  CY * Q * p.S ;
		 -CL * Q * p.S ];

C_bs = [ cos(alpha)   0   -sin(alpha) ;
             0        1        0      ;
     	 sin(alpha)   0    cos(alpha) ];

FA_b = C_bs * FA_s;


%% ---- 4. Aerodynamic moments --------------------------------------------
eta = [ p.Cl_beta * beta ;
		p.Cm0 - p.nt*(p.St*p.lt/(p.S*p.cbar))*(alpha - epsilon) ;
		p.Cn_beta * (1 - p.Cn_beta_alpha*alpha) * beta ];

kq  = p.St*p.lt^2/(p.S*p.cbar^2);      % tail scaling of pitch damping
kdE = p.St*p.lt/(p.S*p.cbar);          % tail scaling of elevator power

dCMdx = (p.cbar/Va) * [ p.Cl_p  p.Cl_q      p.Cl_r ;
                        p.Cm_p  p.Cm_q*kq   p.Cm_r ;
                        p.Cn_p  p.Cn_q      p.Cn_r ];

dCMdu = [ p.Cl_dA  p.Cl_dE      p.Cl_dR ;
          p.Cm_dA  p.Cm_dE*kdE  p.Cm_dR ;
          p.Cn_dA  p.Cn_dE      p.Cn_dR ];

CMac_b = eta + dCMdx*wbe_b + dCMdu*[uA; uE; uR];

MAac_b = CMac_b * Q * p.S * p.cbar;    % dimensional moment [N*m]

rcg_b = [p.Xcg; p.Ycg; p.Zcg];
rac_b = [p.Xac; p.Yac; p.Zac];
MAcg_b = MAac_b + cross(FA_b, rcg_b - rac_b);


%% ---- 5. Engine forces and moments --------------------------------------
F1 = uTH1 * p.m * p.g;      F2 = uTH2 * p.m * p.g;     % [N]
FE1_b = [F1; 0; 0];         FE2_b = [F2; 0; 0];        % thrust along body x
FE_b  = FE1_b + FE2_b;

mew1 = [p.Xcg - p.Xapt1; p.Yapt1 - p.Ycg; p.Zcg - p.Zapt1];
mew2 = [p.Xcg - p.Xapt2; p.Yapt2 - p.Ycg; p.Zcg - p.Zapt2];
MEcg_b = cross(mew1, FE1_b) + cross(mew2, FE2_b);


%% ---- 6. Gravity --------------------------------------------------------
Fg_b = p.m*p.g*[ -sin(theta) ;
				 cos(theta)*sin(phi) ;
				 cos(theta)*cos(phi) ];


%% ---- 7. Equations of motion --------------------------------------------
F_b = Fg_b + FE_b + FA_b;

x1to3_dot = (1/p.m)*F_b - cross(wbe_b, V_b);

x4to6_dot = p.Iinv * (MAcg_b + MEcg_b - cross(wbe_b, p.I*wbe_b));

H_phi = [ 1,  sin(phi)*tan(theta),  cos(phi)*tan(theta) ;
          0,  cos(phi),            -sin(phi)            ;
          0,  sin(phi)/cos(theta),  cos(phi)/cos(theta) ];

x7to9_dot = H_phi * wbe_b;


%% ---- 8. Assemble -------------------------------------------------------
xdot = [x1to3_dot; x4to6_dot; x7to9_dot];

end
