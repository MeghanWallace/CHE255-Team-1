% Peng–Robinson (PR) EOS for methane (CH4) at P = 200 bar
% Units: P in bar, V in L/mol, a in L^2·bar/mol^2, b in L/mol
% R = 0.08314472 L·bar/(mol·K)

clear; clf; close all;

% Constants
P = 200;                % pressure in bar
Tmin = 300; Tmax = 1000; % K
nT = 701;               % number of temperature points
R = 0.08314472;         % L·bar/(mol·K)
Tc = 190.6;             % K
Pc = 45.99;             % bar
omega = 0.011;          % acentric factor for methane

% Parameters
a = 0.45724 * (R^2 * Tc^2) / Pc;   % L^2·bar/mol^2
b = 0.07780 * R * Tc / Pc;         % L/mol

% Temperature vector
T = linspace(Tmin, Tmax, nT);

% Preallocate
V = nan(size(T));
tol_imag = 1e-8;

% For loop
for i = 1:length(T)
    Ti = T(i);

    % Alpha function (Soave-type correction)
    m = 0.37464 + 1.54226*omega - 0.26992*omega^2;
    alpha = (1 + m*(1 - sqrt(Ti/Tc)))^2;

    % Temperature-dependent a(T)
    aT = a * alpha;

    % Cubic Coefficients
    %   P*V^3 - R*T*V^2 + (aT - 3*P*b^2 - 2*R*T*b)*V 
    % + (P*b^3 + R*T*b^2 - aT*b) = 0
    c3 = P;
    c2 = -R * Ti;
    c1 = aT - 3*P*b^2 - 2*R*Ti*b;
    c0 = P*b^3 + R*Ti*b^2 - aT*b;
    coeffs = [c3, c2, c1, c0];

    % Solve cubic for roots
    rootsV = roots(coeffs);

    % Keep only real roots
    real_roots = rootsV(abs(imag(rootsV)) < tol_imag);

    if isempty(real_roots)
        V(i) = NaN;
    else
        V(i) = max(real_roots);
    end
end

% Plot
figure(1);
plot(T, V, 'm-', 'LineWidth', 1.6);
xlabel('Temperature (K)');
ylabel('Molar volume V_m (L/mol)');
title(sprintf('Methane (CH_4) — Peng–Robinson EOS at P = %.0f bar', P));
grid on;
xlim([Tmin Tmax]);

% Sample results
disp('Sample outputs (T, V_m):');
for idx = round(linspace(1,length(T),7))
    fprintf('T = %6.1f K: V_m = %12.6e L/mol\n', T(idx), V(idx));
end
