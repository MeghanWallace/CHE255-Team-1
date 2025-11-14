% Van der Waals EOS for methane (CH4) at 200 bar
% Plots molar volume vs. temperature (300–1000 K)

clear; clf; close all;

% Constants
P = 200;                % bar
Tmin = 300; Tmax = 1000; % K
nT = 701;
R = 0.08314472;         % L·bar/(mol·K)
Tc = 190.6;             % K
Pc = 45.99;             % bar

% Parameters
a = 27 * (R^2 * Tc^2) / (64 * Pc);   % L^2·bar/mol^2
b = R * Tc / (8 * Pc);               % L/mol

% Temperature vector
T = linspace(Tmin, Tmax, nT);
V = nan(size(T));
tol_imag = 1e-8;

% For Loop
for i = 1:length(T)
    Ti = T(i);

    % Cubic Coefficients: P*V^3 - (P*b + R*T)*V^2 + (R*T*b - a)*V + a*b = 0
    c3 = P;
    c2 = -(P*b + R*Ti);
    c1 = R*Ti*b - a;
    c0 = a*b;
    coeffs = [c3, c2, c1, c0];

    rootsV = roots(coeffs);

    % Keep only real & positive roots
    real_roots = real(rootsV(abs(imag(rootsV)) < tol_imag));
    real_roots = real_roots(real_roots > b);

    if isempty(real_roots)
        V(i) = NaN;
    else
        V(i) = max(real_roots);
    end
end

% Plot
figure(1);
plot(T, V, 'b-', 'LineWidth', 1.6);
xlabel('Temperature (K)');
ylabel('Molar volume V_m (L/mol)');
title(sprintf('Methane (CH_4) — Van der Waals EOS at P = %.0f bar', P));
grid on;
xlim([Tmin Tmax]);
ylim([0.1,0.5]);

% Sample results
disp('Sample outputs (T, V_m):');
for idx = round(linspace(1,length(T),7))
    fprintf('T = %6.1f K: V_m = %12.6e L/mol\n', T(idx), V(idx));
end
