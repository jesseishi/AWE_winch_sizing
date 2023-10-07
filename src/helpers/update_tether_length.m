function [E_eff, CR_eff, C] = update_tether_length(tether_length, kite, tether, environment)
% calculate the effective lift-to-drag ratio of a kite, taking
% the tether drag into account.
CD = kite.CL / kite.E;
CD_eff = CD + 0.5 * tether.r_m * tether_length / kite.S_m2 * tether.Cd_t;
E_eff = kite.CL / CD_eff;

% Calculate the effective resultant coefficient and C (mathscr{C}) related
% to the performance of the kite.
CR_eff = kite.CL * sqrt(1 + 1/E_eff^2);
C = 0.5 * environment.rho_kgpm3 * kite.S_m2 * CR_eff * (1 + E_eff^2);
end
