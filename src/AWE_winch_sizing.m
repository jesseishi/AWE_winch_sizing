%% Load parameters.
close all; clear; clc;

load MegAWES
addpath helpers

Lt_m = 1000;  % tether length.
[kite.E_eff, kite.CR_eff, kite.C] = update_tether_length(Lt_m, kite, tether, environment);
winch.K_w = winch.J_kgm2 / winch.r_m^2;

%% Trim analysis.
% Where f = f^* so F_t = tau * r.
figure
v_r = 0:0.1:7;
plot(v_r, 4 * kite.C * v_r.^2 ./ 1e6, 'DisplayName', 'winch control curve')
hold on
xlabel('Reel-out speed (m/s)')
ylabel('Tether force (MN)')
grid on
xlim([0, max(v_r)])

saveas(gcf, '../Results/vrFt', 'epsc')

%% Build the system.
% Set trim condition.
v_w0 = 10;  % wind speed, m/s (achieves a nice scale on the winch control curve later)

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

%% Normalize the transfer functions.
vrn_over_vw = vr_over_vw / dcgain(vr_over_vw);
fn_over_vw = f_over_vw;  % dcgain = 0.
Ftn_over_vw = 36 / (32 * kite.C * v_w0) * Ft_over_vw;
Ftn_over_vw2 = Ft_over_vw / dcgain(Ft_over_vw);  % correct
Pn_over_vw = P_over_vw / dcgain(P_over_vw);

%% Bode plot.
figure
opts = bodeoptions;
opts.MagUnits = 'abs';
bodeplot(vrn_over_vw, Ftn_over_vw, opts)
legend('reel-out speed (normalized)', 'tether force (normalized)')
grid
print('../Results/winch_bode.eps','-depsc2');

%% Show response on winch control curve.
% Plot the winch control curve.
figure
v_r0 = v_w0/3;
v_r = 0:0.1:v_r0*1.4;
plot(v_r, 4 * kite.C * v_r.^2 / 1e6, 'k--', 'DisplayName', 'winch control curve')
hold on
xlabel('Reel-out speed (m/s)')
ylabel('Tether force (MN)')
set(gca,'ColorOrderIndex',1)


F_t0 = 4/9 * kite.C * v_w0^2 / 1e6;

% For different input frequencies, plot the response.
w_0s = [1e2, 1e3, 1e4];
H = [vr_over_vw; Ft_over_vw];
for i = 1:length(w_0s)
    w_0 = w_0s(i);  % input frequency.
    T = 2*pi/w_0;  % period.
    
    % Run the simulation
    t = 0:T/100:3*T;
    u = 2*sin(w_0*t);
    y = lsim(H, u, t);
    
    % Plot the response.
    plot(v_r0 + y(:, 1), F_t0 + y(:, 2)./1e6, 'DisplayName', sprintf('w_0 = 10^%d rad/s', log10(w_0)))
end
legend('Location','NorthWest')
grid on

print('../Results/linear_response.eps','-depsc2');

%% Absolute value of the normalized tether force overshoot.
% make a plot with trim apparent wind speed on the x-axis, K_w on the
% y-azis and the normalized tether force in color.
w_0 = pi/10;  % Apparent wind speed oscillation of MegAWES.
w_0 = w_0*3;

N = 67;
K_wv = logspace(3.5, 6.5, N);
v_w0v = logspace(log10(0.01), log10(30), N);
[K_W, v_w0] = meshgrid(K_wv, v_w0v);

MAGS = sqrt((9 * K_W.^2 .* w_0^2 + 64 * kite.C^2 * v_w0.^2) ./ ...
    (4 * K_W.^2 * w_0^2 + 64 * kite.C^2 * v_w0.^2));


figure
contourf(v_w0, K_W, MAGS, linspace(1, 1.5, 101), 'EdgeColor', 'none');
hold on
contour(v_w0, K_W, MAGS, [1.01, 1.05, 1.10], '-k', 'ShowText', true, "LabelFormat", "%.2f")

colormap(flipud(parula))
clim([1.0, 1.5])  % This can really mess with the plot. Double check that it does not only change the scale of the legend and not on the plot!
c = colorbar;
c.Label.String = 'Normalized tether force (-)';
xlabel('Trim equivalent wind speed (m/s)')
ylabel('Winch sizing parameter (kg)')
% set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

% print('../Results/MegAWES_sizing.eps','-depsc2');
set(gcf,'renderer','opengl');   % For avoiding vertical/horizontal line artefacts when saving to eps.
saveas(gcf, '../Results/MegAWES_sizing', 'epsc')

%% Closed-form solution
w_0 = pi/10;
w_0 = w_0 * 3;
v_w0 = 10;  % cut-in wind speed
Ft_overshoot = 1.01;  % Requirement

K_w_select = sqrt((64 * kite.C^2 * v_w0^2 * (Ft_overshoot^2 - 1)) / ...
    ((9 - 4 * Ft_overshoot^2) * w_0^2));

r_select = 2.0;
J_select = K_w_select * r_select^2;



