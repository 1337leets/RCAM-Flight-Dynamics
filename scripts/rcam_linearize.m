function [A, B, Alat, Blat, Alon, Blon] = rcam_linearize(x_tr, u_tr, p)
%RCAM_LINEARIZE  Numerical linearization about a trim point.
%
%   [A, B, Alat, Blat, Alon, Blon] = rcam_linearize(x_tr, u_tr, p)
%
%   xdot = f(x,u) is nonlinear, but close to the trim point any smooth
%   function looks linear (first-order Taylor expansion):
%
%       xdot ~ f(x0,u0) + (df/dx)|_0 * dx + (df/du)|_0 * du
%              \_______/    \______/          \______/
%               = 0 (trim)      A                B
%
%   A and B are the Jacobians of f. Instead of deriving them by hand, they
%   are computed numerically with central differences:
%
%       df/dx_j ~ [ f(x + h*e_j) - f(x - h*e_j) ] / (2h)
%
%   In other words, each state is nudged in turn and the change in xdot is
%   recorded. Because this works at any flight condition, it is also the
%   basis for gain scheduling.

if nargin < 3, p = rcam_params(); end

nx = numel(x_tr);
nu = numel(u_tr);

A = zeros(nx, nx);
B = zeros(nx, nu);

% Step size. Too small amplifies round-off noise, too large lets the
% nonlinearity leak in. A fixed step of 1e-6 works well for this model.
h = 1e-6;

for j = 1:nx
    dx = zeros(nx,1);  dx(j) = h;
    f_plus  = rcam_dynamics(x_tr + dx, u_tr, p);
    f_minus = rcam_dynamics(x_tr - dx, u_tr, p);
    A(:,j) = (f_plus - f_minus) / (2*h);
end

for j = 1:nu
    du = zeros(nu,1);  du(j) = h;
    f_plus  = rcam_dynamics(x_tr, u_tr + du, p);
    f_minus = rcam_dynamics(x_tr, u_tr - du, p);
    B(:,j) = (f_plus - f_minus) / (2*h);
end

% ---- Extract the subsystems ---------------------------------------------
% In symmetric flight (beta = 0, phi = 0) the longitudinal and
% lateral-directional motions decouple, so the 9x9 problem splits into
% two independent smaller ones.
%
% Lateral-directional states:   v(2)  p(4)  r(6)  phi(7)
% Lateral-directional controls: uA(1)  uR(3)
ilat = [2 4 6 7];
jlat = [1 3];
Alat = A(ilat, ilat);
Blat = B(ilat, jlat);

% Longitudinal states:   u(1)  w(3)  q(5)  theta(8)
% Longitudinal controls: uE(2)  uTH1(4)
% Note: Blon holds the left engine only. For symmetric thrust use
% B(ilon,4) + B(ilon,5).
ilon = [1 3 5 8];
jlon = [2 4];
Alon = A(ilon, ilon);
Blon = B(ilon, jlon);

end
