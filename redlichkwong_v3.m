% Redlich Kwong EOS for methane (CH4) at 200 bar
% Units: P in bar, V in L/mol, a in L^2·bar·K^0.5/mol^2, b in L/mol
% R = 0.08314472 L·bar/(mol·K)

clear; clf; close all;

% Constants
P = 200;                % pressure in bar
Tmin = 300; Tmax = 1000; % K
nT = 701;               % number of T points
R = 0.08314472;         % L·bar/(mol·K)
Tc = 190.6;             % K
Pc = 45.99;             % bar

% Redlich-Kwong parameters (from Tc, Pc)
% a has units L^2·bar·K^0.5/mol^2, b has units L/mol
a = 0.42748 * R^2 * Tc^(2.5) / Pc;
b = 0.08664 * R * Tc / Pc;

% Temperature vector
T = linspace(Tmin, Tmax, nT);

% Preallocate
V = nan(size(T));        
tol_imag = 1e-8;          % tolerance for imaginary parts

% For loop
for i = 1:length(T)
    Ti = T(i);
    A = a / sqrt(Ti);   % RK temperature-dependent term

    % Cubic coefficients: P*V^3 - R*T*V^2 + (-P*b^2 - R*T*b + A)*V - A*b = 0
    c3 = P;
    c2 = -R * Ti;
    c1 = -P * b^2 - R * Ti * b + A;
    c0 = -A * b;
    coeffs = [c3, c2, c1, c0];

    rootsV = roots(coeffs);

    % keep only real roots
    real_roots = rootsV(abs(imag(rootsV)) < tol_imag);
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
title(sprintf('Methane (CH_4) — Redlich–Kwong EOS at P = %.0f bar', P));
grid on;
xlim([Tmin Tmax]);

% Sample outputs
disp('Sample outputs (T, V_m):');
for idx = round(linspace(1,length(T),7))
    fprintf('T = %6.1f K: V_m = %12.6e L/mol\n', T(idx), V(idx));
end
