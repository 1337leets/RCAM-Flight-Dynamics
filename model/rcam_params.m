function p = rcam_params()
%RCAM_PARAMS  GARTEUR Research Civil Aircraft Model parameter set.
%
%   Source: GARTEUR FM(AG-08), "Robust Flight Control Design Challenge
%   Problem Formulation and Manual: the Research Civil Aircraft Model
%   (RCAM)", TP-088-3, 1995. The values follow the February 1997
%   revision as used in common open implementations. The moment
%   coefficients were checked against TP-088-3 Version 1, equation 2.33
%   (page 16).
%
%   Aircraft: a twin-engine wide-body civil transport of about 120 t
%   (A310/B767 class).

% ---- Mass and geometry --------------------------------------------------
p.m      = 120000;      % kg    aircraft mass
p.cbar   = 6.6;         % m     mean aerodynamic chord
p.S      = 260;         % m^2   wing area
p.St     = 64;          % m^2   horizontal tail area
p.lt     = 24.8;        % m     tail moment arm

% ---- Reference points (normalized by cbar) ------------------------------
p.Xcg = 0.23*p.cbar;  p.Ycg = 0.0;  p.Zcg = 0.10*p.cbar;   % centre of gravity
p.Xac = 0.12*p.cbar;  p.Yac = 0.0;  p.Zac = 0.0;           % aerodynamic centre

% ---- Engine locations (body axes, m) ------------------------------------
p.Xapt1 = 0;  p.Yapt1 = -7.94;  p.Zapt1 = -1.9;   % left engine
p.Xapt2 = 0;  p.Yapt2 =  7.94;  p.Zapt2 = -1.9;   % right engine

% ---- Inertia tensor (per unit mass, m^2) --------------------------------
%   Actual inertia: p.I = p.m * p.Ib
p.Ib = [ 40.07,  0,      -2.0923;
          0,     64.0,    0;
         -2.0923, 0,      99.92 ];
p.I    = p.m * p.Ib;
p.Iinv = inv(p.I);

% ---- Environment --------------------------------------------------------
p.rho = 1.225;          % kg/m^3  (sea level; override for other altitudes)
p.g   = 9.81;           % m/s^2

% ---- Lift ---------------------------------------------------------------
p.depsda      = 0.25;                 % downwash gradient d(eps)/d(alpha)
p.alpha_L0    = -11.5*pi/180;         % rad    zero-lift angle of attack
p.alpha_switch= 14.5*pi/180;          % rad    switch from linear to stall model
p.n           = 5.5;                  % 1/rad  wing lift slope (linear region)
% Cubic fit beyond the switch point: CL_wb = a3*a^3 + a2*a^2 + a1*a + a0
p.a3 = -768.5;  p.a2 = 609.2;  p.a1 = -155.2;  p.a0 = 15.212;
p.nt = 3.1;                           % 1/rad  tail lift slope

% ---- Drag ---------------------------------------------------------------
p.CDmin = 0.13;
p.d1    = 0.07;
p.d0    = 0.654;
% Form used in the RCAM manual:  CD = CDmin + d1*(n*alpha + d0)^2
%   so drag grows with the square of the angle of attack (induced drag).

% ---- Side force ---------------------------------------------------------
p.CY_beta = -1.6;       % 1/rad
p.CY_dR   =  0.24;      % 1/rad   rudder to side force

% ---- Static moment coefficients -----------------------------------------
p.Cl_beta = -1.4;       % 1/rad   dihedral effect (negative is stable)
p.Cm0     = -0.59;
p.Cm_alpha= -3.1;       % 1/rad   longitudinal static stability (negative)
p.Cn_beta =  1.0;       % 1/rad   weathercock stability at alpha = 0 (positive)
p.Cn_beta_alpha = 180/(15*pi);  % 1/rad  Cn = Cn_beta*(1 - Cn_beta_alpha*alpha)*beta
                        %   (directional stability vanishes at alpha = 15 deg)

% ---- Rate derivatives ---------------------------------------------------
p.Cl_p = -11.0;   p.Cl_q =  0.0;    p.Cl_r =  5.0;
p.Cm_p =   0.0;   p.Cm_q = -4.03;   p.Cm_r =  0.0;
p.Cn_p =   1.7;   p.Cn_q =  0.0;    p.Cn_r = -11.5;

%   Note: Cm_q and Cm_dE are normalized to the tail. rcam_dynamics scales
%   them by St*lt^2/(S*cbar^2) and St*lt/(S*cbar), as in TP-088-3 eq. 2.33.

% ---- Control derivatives ------------------------------------------------
p.Cl_dA = -0.6;   p.Cl_dE =  0.0;   p.Cl_dR =  0.22;
p.Cm_dA =  0.0;   p.Cm_dE = -3.1;   p.Cm_dR =  0.0;
p.Cn_dA =  0.0;   p.Cn_dE =  0.0;   p.Cn_dR = -0.63;

% ---- Control limits (rad, RCAM manual) ----------------------------------
p.uA_lim = [-25, 25] * pi/180;      % aileron
p.uE_lim = [-25, 10] * pi/180;      % elevator (all-moving tailplane)
p.uR_lim = [-30, 30] * pi/180;      % rudder
p.uTH_lim= [ 0.5, 10] * pi/180;     % throttle (an angle in RCAM)

% ---- Actuator time constants (s) ----------------------------------------
p.tau_A = 0.15;
p.tau_E = 0.15;
p.tau_R = 0.30;         % The rudder is slower than the other surfaces,
                        % which limits how much a yaw damper can achieve.

end
