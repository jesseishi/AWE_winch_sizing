% MegAWES parameters obtained from eijkelhofSixdegreesoffreedomSimulationModel2022.

kite = struct( ...
    'CL', 1.8, ...           % lift coefficient [-]
    'E', 8.4, ...            % lift-drag ratio exclusive tether [-], qsm.py will calculate tether contribution.
    'S_m2', 150.45, ...      % wing reference area [m^2]
    'm_kg', 6885.2);         % mass [kg]

tether = struct( ...
    'rho_kgpm3', 242.8, ...  % tether density [kg/m3]
    'r_m', 0.01485, ...      % tether radius [m]
    'Cd_t', 1.1);            % tether drag coefficient [-]

winch = struct( ...
    'r_m', 0.4, ...          % radius [m]
    'J_kgm2', 32.0, ...      % inertia [kgm^2]
    'friction', 10.0);       % rotational friction coefficient [Nm/(rad/s)]

environment = struct( ...
    'rho_kgpm3', 1.225);     % density [kg/m^3]

save("MegAWES", "kite", "tether", "winch", "environment");
