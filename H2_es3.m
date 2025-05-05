close all
clear all
a=1;
r=1;
w_rand=rand(1,1)*pi*2;

out = sim('HW2_es3');


% Impostazioni generali
figure('Color', 'w'); % sfondo bianco

% === GRAFICO 1: Posizione punto B e centro della ruota ===
subplot(3,1,1)
plot(out.out1(:,4), out.out1(:,5), 'b', 'LineWidth', 1.5); hold on;
plot(out.out1(:,9), out.out1(:,10), 'g', 'LineWidth', 1.5);
plot(out.out1(:,6), out.out1(:,7), 'r--', 'LineWidth', 1.5);

legend('linear control trajectory', 'almost non-linear control trajectory', 'desired trajectory');
title('Position in the plan');
xlabel('X [m]');
ylabel('Y [m]');
grid on;

% === GRAFICO 2: Errori ===
subplot(3,1,2)
plot(out.out1(:,8), out.out1(:,1), 'g', 'LineWidth', 1.5); hold on;
plot(out.out1(:,8), out.out1(:,2), 'm', 'LineWidth', 1.5);
plot(out.out1(:,8), out.out1(:,3), 'b', 'LineWidth', 1.5);
legend('error 1', 'error 2','error theta');
title('tracking errors linear control ');
xlabel('Time [s]');
ylabel('error');
grid on;

% === GRAFICO 3: Errori ===
subplot(3,1,3)
plot(out.out1(:,8), out.out1(:,11), 'g', 'LineWidth', 1.5); hold on;
plot(out.out1(:,8), out.out1(:,12), 'm', 'LineWidth', 1.5);
plot(out.out1(:,8), out.out1(:,13), 'b', 'LineWidth', 1.5);
legend('error 1', 'error 2','error theta');
title('tracking errors Almost non-linear control ');
xlabel('Time [s]');
ylabel('error');
grid on;

distanza1 = sqrt( (out.out1(round(length(out.out1(:,8))/10):length(out.out1(:,8)),4) - out.out1(round(length(out.out1(:,8))/10):length(out.out1(:,8)),6)).^2 + (out.out1(round(length(out.out1(:,8))/10):length(out.out1(:,8)),5) - out.out1(round(length(out.out1(:,8))/10):length(out.out1(:,8)),7)).^2 );

% Traiettoria 2 vs Riferimento
distanza2 = sqrt( (out.out1(round(length(out.out1(:,8))/10):length(out.out1(:,8)),9) - out.out1(round(length(out.out1(:,8))/10):length(out.out1(:,8)),6)).^2 + (out.out1(round(length(out.out1(:,8))/10):length(out.out1(:,8)),10) - out.out1(round(length(out.out1(:,8))/10):length(out.out1(:,8)),7)).^2 );



% Distanza massima
max_distanza1 = max(distanza1);
max_distanza2 = max(distanza2);

% Output a video
fprintf('Maximum distance between Trajectory linear control and Reference: %.4f m\n', max_distanza1);
fprintf('Maximum distance between Trajectory Almost non-linear and Reference: %.4f m\n', max_distanza2);