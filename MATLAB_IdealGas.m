% Ideal Gas Law for methane (CH4) at 200 bar
% constants
P = 200;                   % bar
R = 0.08314;               % L·bar/(mol·K) 
Tmin = 300; Tmax = 1000;   % K

% ideal gas law
T = 300:50:1000;
V = (R*T)/P;
plot(T,V, 'LineWidth', 1.6, 'DisplayName', 'Ideal Gas');


% plot apperance
figure(1);
hold on;

xlabel('Temperature (K)');
ylabel('Molar Volume V_m (L/mol)');
title(sprintf('Methane: Ideal Gas V_m vs T at P = %.0f bar', P));
grid on;
xlim([Tmin Tmax]);
ylim([0.1,0.5]);
hold off;

% Sample results
disp('Sample outputs (T, V_m):');
for idx = round(linspace(1,length(T),7))
    fprintf('T = %6.1f K: V_m = %12.6e L/mol\n', T(idx), V(idx));
end
