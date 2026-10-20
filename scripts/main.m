%MAIN  RCAM driver script: trim, linearize and list the flight modes.
%
%  Run from the scripts/ folder, or add the repository to the path first
%  with addpath(genpath(<repo root>)).

clear; clc;
addpath('../model'); addpath('.'); addpath('../tests');

p  = rcam_params();
Va = 85;                       % m/s  (about 165 kt)

[x_tr, u_tr, info] = rcam_trim(Va, p);

[A, B, Alat, Blat, Alon, Blon] = rcam_linearize(x_tr, u_tr, p);

fprintf('\nLateral-directional A matrix (states: v, p, r, phi):\n');
disp(Alat);
fprintf('Lateral-directional B matrix (controls: aileron, rudder):\n');
disp(Blat);

analyze_modes(Alat, Alon);
