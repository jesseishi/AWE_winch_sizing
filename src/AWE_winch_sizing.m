%% Load parameters.
close all; clear; clc;

load MegAWES
addpath helpers

Lt_m = 1000;  % tether length.
[kite.E_eff, kite.CR_eff, kite.C] = update_tether_length(Lt_m, kite, tether, environment);
winch.K_w = winch.J_kgm2 / winch.r_m^2;

%% Trim analysis.
% Where f = f^* so F_t = tau.
% TODO.

%% Build the system.
% Set trim condition.
v_w0 = 5;  % wind speed, m/s

vr_over_vw = tf(4 * kite.C * v_w0, [3 * winch.K_w, 12 * kite.C * v_w0]);

% Derived on paper (see notebook).
f_over_vw = tf([-winch.K_w, 0], [3 * v_w0 * winch.K_w, 12 * kite.C * v_w0^2]);
Ft_over_vw = tf([12 * kite.C * v_w0 * winch.K_w, 32 * kite.C^2 * v_w0^2], ...
    [9 * winch.K_w, 36 * kite.C * v_w0]);
P_over_vw = tf(16 * kite.C^2 * v_w0^3, [9 * winch.K_w, 36 * kite.C * v_w0]);

% Derived half symbolically and half analytically.
f_over_vw2 = 1/v_w0 * vr_over_vw - 1/(3*v_w0);  % correct
Ft_over_vw2 = 4/3 * kite.C * v_w0 * (1 - vr_over_vw);  % correct
P_over_vw2 = 4/3 * kite.C * v_w0^2 * vr_over_vw;  % correct

%% Analysis.
% Normalize the transfer functions.
vrn_over_vw = vr_over_vw / dcgain(vr_over_vw);
fn_over_vw = f_over_vw;  % dcgain is 0.
Ftn_over_vw = Ft_over_vw / dcgain(Ft_over_vw);
Pn_over_vw = P_over_vw / dcgain(P_over_vw);

opts = bodeoptions;
opts.MagUnits = 'abs';
opts.PhaseVisible = 'off';

bodeplot(vrn_over_vw, fn_over_vw, Ftn_over_vw, Pn_over_vw, opts)
legend
grid
