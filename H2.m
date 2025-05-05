close all
clear all

b=1;
%inizializzazione
qi =[0, 0, 0];
xi = qi(1);
yi = qi(2);
theta_i=qi(3);
% % 
% qf = rand(1,3)';
% qf = qf / norm(qf); 

 qf = [0.7004
    0.1144
    0.7045];
xf = qf(1);
yf = qf(2);
theta_f = qf(3);

k = 1;
T=1;
t = 0:0.01:T;

    a_2 = 3 / T^2;
    a_3 = -2 / T^3;

    % Calcolo della traiettoria 
    s = a_3 * t.^3 + a_2 * t.^2;
    ds = 3 * a_3 * t.^2 + 2 * a_2 * t; 
    dds = 6 * a_3 * t + 2 * a_2; 
  
       
   
    

    alpha_x = k*cos(theta_f)-3*xf; 
    beta_x = k*cos(theta_i)+3*xi;
    x_s = s.^3.*qf(1) - (s-1).^3.*qi(1) + alpha_x.*s.^2.*(s-1) + beta_x.*s.*(s-1).^2;
    dx = 3 .* s.^2 .* qf(1) - 3 .* (s - 1).^2 .* qi(1) + alpha_x .* s .* (3 .* s - 2) + beta_x .* (3 .* s - 1) .* (s - 1);
    ddx = 6 .* s .* qf(1) - 6 .* (s - 1) .* qi(1) + alpha_x .* (6 .* s - 2) + beta_x .* (6 .* s - 4);

    alpha_y = k*sin(theta_f)-3*yf;
    beta_y=  k*sin(theta_i)+3*yi;
    y_s = s.^3.*qf(2) - (s-1).^3.*qi(2) + alpha_y.*s.^2.*(s-1) + beta_y.*s.*(s-1).^2;
    dy = 3 .* s.^2 .* qf(2) - 3 .* (s - 1).^2 .* qi(2) + alpha_y .* s .* (3 .* s - 2) + beta_y .* (3 .* s - 1) .* (s - 1);
    ddy = 6 .* s .* qf(2) - 6 .* (s - 1) .* qi(2) + alpha_y .* (6 .* s - 2) + beta_y .* (6 .* s - 4);

    theta = atan2(dy, dx);
    v_tilde = sqrt(dx.^2 + dy.^2);
    w_tilde = (ddy .* dx - ddx .* dy) ./ ((dx).^2 + (dy).^2);

    % Calcolo di v_t e w_t
    v_t = v_tilde .* ds;
    w_t = w_tilde .* ds;


%% =========================================================================
%  PLOT 1: TRAIETTORIA SPAZIALE (x vs y)
%% =========================================================================


subplot(3,2,1);
plot(x_s, y_s, 'b-', 'LineWidth', 2); hold on;
plot(0, 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % punto iniziale
plot(qf(1), qf(2), 'gx', 'MarkerSize', 8, 'LineWidth', 2); % punto finale
xlabel('x(t)');
ylabel('y(t)');
title('Traiettoria spaziale: x(t) vs y(t)');
legend({'Traiettoria', 'Start', 'Goal'}, 'Location', 'best');
grid on;
axis equal;

%% =========================================================================
%  PLOT 2: s(t), s'(t), s''(t)
%% =========================================================================

subplot(3,2,3);
plot(t, s,    'k',  'LineWidth', 2); hold on;
plot(t, ds,   'k--', 'LineWidth', 2);
plot(t, dds,  'k-.', 'LineWidth', 2);
xlabel('t');
ylabel('Valori');
title('s(t), s''(t), s''''(t) nel tempo');
legend({'s(t)', 's''(t)', 's''''(t)'}, 'Location', 'best');
grid on;

%% =========================================================================
%  PLOT 3: v(t) e w(t)
%% =========================================================================

subplot(3,2,5);
plot(t, v_t, 'g-', 'LineWidth', 2); hold on;
plot(t, w_t, 'r-', 'LineWidth', 2);
xlabel('t');
ylabel('Velocità');
title('Velocità lineare v(t) e angolare w(t)');
legend({'v(t)', 'w(t)'}, 'Location', 'best');
grid on;

v_max = 0.5; 
w_max = 2; 

scaling_value=max(max(abs(v_t)/v_max),max(abs(w_t))/w_max);
T=T*scaling_value;
a2 = 3/T^2;
a3 = -2/T^3;

t = 0:0.01:T;

    a_2 = 3 / T^2;
    a_3 = -2 / T^3;

    % Calcolo della traiettoria 
    s = a_3 * t.^3 + a_2 * t.^2;
    ds = 3 * a_3 * t.^2 + 2 * a_2 * t; 
    dds = 6 * a_3 * t + 2 * a_2; 
  
    alpha_x = k*cos(theta_f)-3*xf;
    beta_x = k*cos(theta_i)+3*xi;
    x_s = s.^3.*qf(1) - (s-1).^3.*qi(1) + alpha_x.*s.^2.*(s-1) + beta_x.*s.*(s-1).^2;
    dx = 3 .* s.^2 .* qf(1) - 3 .* (s - 1).^2 .* qi(1) + alpha_x .* s .* (3 .* s - 2) + beta_x .* (3 .* s - 1) .* (s - 1);
    ddx = 6 .* s .* qf(1) - 6 .* (s - 1) .* qi(1) + alpha_x .* (6 .* s - 2) + beta_x .* (6 .* s - 4);

    alpha_y = k*sin(theta_f)-3*yf;
    beta_y=  k*sin(theta_i)+3*yi;
    y_s = s.^3.*qf(2) - (s-1).^3.*qi(2) + alpha_y.*s.^2.*(s-1) + beta_y.*s.*(s-1).^2;
    dy = 3 .* s.^2 .* qf(2) - 3 .* (s - 1).^2 .* qi(2) + alpha_y .* s .* (3 .* s - 2) + beta_y .* (3 .* s - 1) .* (s - 1);
    ddy = 6 .* s .* qf(2) - 6 .* (s - 1) .* qi(2) + alpha_y .* (6 .* s - 2) + beta_y .* (6 .* s - 4);

    theta = atan2(dy, dx);
    v_tilde = sqrt(dx.^2 + dy.^2);
    w_tilde = (ddy .* dx - ddx .* dy) ./ ((dx).^2 + (dy).^2);

    % Calcolo di v_t e w_t
    v_t = v_tilde .* ds;
    w_t = w_tilde .* ds;






%% =========================================================================
%  PLOT 1: TRAIETTORIA SPAZIALE (x vs y)
%% =========================================================================


subplot(3,2,2);
plot(x_s, y_s, 'b-', 'LineWidth', 2); hold on;
plot(0, 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % punto iniziale
plot(qf(1), qf(2), 'gx', 'MarkerSize', 8, 'LineWidth', 2); % punto finale
xlabel('x(t)');
ylabel('y(t)');
title('Traiettoria spaziale SCALATA: x(t) vs y(t)');
legend({'Traiettoria', 'Start', 'Goal'}, 'Location', 'best');
grid on;
axis equal;

%% =========================================================================
%  PLOT 2: s(t), s'(t), s''(t)
%% =========================================================================

subplot(3,2,4);
plot(t, s,    'k',  'LineWidth', 2); hold on;
plot(t, ds,   'k--', 'LineWidth', 2);
plot(t, dds,  'k-.', 'LineWidth', 2);
xlabel('t');
ylabel('Valori');
title('s(t), s''(t), s''''(t)SCALATA nel tempo');
legend({'s(t)', 's''(t)', 's''''(t)'}, 'Location', 'best');
grid on;

%% =========================================================================
%  PLOT 3: v(t) e w(t)
%% =========================================================================

subplot(3,2,6);
plot(t, v_t, 'g-', 'LineWidth', 2); hold on;
plot(t, w_t, 'r-', 'LineWidth', 2);
xlabel('t');
ylabel('Velocità');
title('Velocità lineare v(t) e angolare w(t) SCALATA');
legend({'v(t)', 'w(t)'}, 'Location', 'best');
grid on;






%% ex 2

k1=5;
k2=5;

y1=x_s+b*cos(theta);
y2=y_s+b*sin(theta);
y1_dot=cos(theta).*v_t-b*sin(theta).*w_t;
y2_dot=sin(theta).*v_t+b*cos(theta).*w_t;


y1_d = timeseries(y1,t);
y2_d = timeseries(y2,t); 
y1_dot_d = timeseries(y1_dot,t); 
y2_dot_d = timeseries(y2_dot,t); 
theta_d = timeseries(theta,t) ;
x_d = timeseries(x_s,t); 
y_d = timeseries(y_s,t); 
theta_t = timeseries(theta,t);

Tsca=T;
x_traj=x_d;
y_traj=y_d;

out = sim('HW2_es2');


% Impostazioni generali
figure('Color', 'w'); % sfondo bianco

% === GRAFICO 1: Posizione punto B e centro della ruota ===
subplot(3,1,1)
plot(out.out1(:,1), out.out1(:,2), 'b', 'LineWidth', 1.5); hold on;
plot(out.out1(:,3), out.out1(:,4), 'r', 'LineWidth', 1.5);
plot(y1,y2, 'b--', 'LineWidth', 1.5);
legend('Point B', 'Wheel Center');
title('trajectory');
xlabel('X [m]');
ylabel('Y [m]');
grid on;

% === GRAFICO 2: Errori ===
subplot(3,1,2)
plot(out.out1(:,9), out.out1(:,5), 'g', 'LineWidth', 1.5); hold on;
plot(out.out1(:,9), out.out1(:,6), 'm', 'LineWidth', 1.5);
plot(out.out1(:,9), out.out1(:,10), 'b', 'LineWidth', 1.5);
legend('Error x', 'Error y','Error theta');
title('Tracking errors');
xlabel('Time [s]');
ylabel('Errors');
grid on;

% === GRAFICO 3: Comandi v e w ===
subplot(3,1,3)
plot(out.out1(:,9), out.out1(:,7), 'c', 'LineWidth', 1.5); hold on;
plot(out.out1(:,9), out.out1(:,8), 'k', 'LineWidth', 1.5);
legend('Linear velocity v', 'Angular velocity w');
title('real values ​​of v e w ');
xlabel('Time [s]');
ylabel('velocity');
grid on;





error_norm = sqrt(out.out1(:,5).^2 + out.out1(:,6).^2 + out.out1(:,10).^2);


plot(out.out1(:,9), error_norm, 'r', 'LineWidth', 2);
title('Norma dell''errore totale');
xlabel('Time [s]');
ylabel('||Error||');
grid on;

hold off
a = load('ERROR_NORM_0.01.mat');
b = load('ERROR_NORM_0.1.mat');
bc= load('ERROR_NORM_0.5.mat');
c = load('ERROR_NORM_1.mat');
d = load('ERROR_NORM_5.mat');
e = load('ERROR_NORM_10.mat');

figure; hold on; grid on;

plot(out.out1(:,9),a.error_norm, 'LineWidth', 1.5);
plot(out.out1(:,9),b.error_norm, 'LineWidth', 1.5);
plot(out.out1(:,9),bc.error_norm, 'LineWidth', 1.5);
plot(out.out1(:,9),c.error_norm, 'LineWidth', 1.5);
plot(out.out1(:,9),d.error_norm, 'LineWidth', 1.5);
plot(out.out1(:,9),e.error_norm, 'LineWidth', 1.5);

legend('b=0.01', 'b=0.1','b=0.5', 'b=1', 'b=5', 'b=10');
title('Comparison Norm Error at Variation of b');
xlabel('Sample');
ylabel('||Error||');


% saveas(gcf, 'Hw_es2_b100cm.png')
% saveas(gcf, 'Hw_es2_b100cm.pdf')