% Ideal Gas Law for methane (CH4) at 200 bar
% constants
P = 200;                   % bar
R = 0.08314;               % L·bar/(mol·K) 
Tmin = 300; Tmax = 1000;   % K

% ideal gas law
T = 300:50:1000;
V1 = (R*T)/P;
plot(T,V1, 'LineWidth', 1.6, 'DisplayName', 'Ideal Gas');


% plot apperance
figure(1);
hold on;

xlabel('Temperature (K)');
ylabel('Molar Volume V_m (L/mol)');
title(sprintf('Methane: Ideal Gas V_m vs T at P = %.0f bar', P));
legend('Location', 'best');
grid on;
xlim([Tmin Tmax]);
hold off
