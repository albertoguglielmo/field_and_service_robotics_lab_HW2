close all
clear all

Ts=0.02;
a=1;
qi=[1+a,1,pi/4];

out = sim('HW2_es4');


% Impostazioni generali
figure('Color', 'w'); % sfondo bianco

% === GRAFICO 1: Posizione punto B e centro della ruota ===
subplot(3,1,1)
plot(out.out1(:,1), out.out1(:,2), 'b', 'LineWidth', 2); hold on;
plot(qi(1), qi(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2); % punto iniziale
plot(0, 0, 'gx', 'MarkerSize', 8, 'LineWidth', 2); % punto finale
xlabel('x(t)');
ylabel('y(t)');
title('Traiettoria spaziale: x(t) vs y(t)');
legend({'Traiettoria', 'Start', 'Goal'}, 'Location', 'best');
grid on;
axis equal;

% === GRAFICO 2: Errori ===
subplot(3,1,2)
plot(out.out1(:,7), out.out1(:,1), 'g', 'LineWidth', 1.5); hold on;
plot(out.out1(:,7), out.out1(:,2), 'm', 'LineWidth', 1.5);
plot(out.out1(:,7), out.out1(:,3), 'b', 'LineWidth', 1.5);
legend('x ', 'y ','theta');
title(' STATES ');
xlabel('Time [s]');
ylabel('state');
grid on;

% === GRAFICO 3: Errori ===
subplot(3,1,3)
stairs(out.out1(:,7), out.out1(:,4), 'g', 'LineWidth', 1.5); hold on;
stairs(out.out1(:,7), out.out1(:,5), 'm', 'LineWidth', 1.5);
stairs(out.out1(:,7), out.out1(:,6), 'b', 'LineWidth', 1.5);
legend('x discretised', 'y discretised','theta discretised');
title(' Runge-Kutta ');
xlabel('Time [s]');
ylabel('state');
grid on;


%saveas(gcf, 'Hw_es4.png')