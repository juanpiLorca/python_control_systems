% Bode:
% Gráficos de acuerdo al modelo corregido por el profesor: Miguel Torres

% posición carro: x
% Valores aproximados para los coeficientes de p(s) en X(s):
num_x = 10^-5*[2.3121 48.1696 -2.8353];
% Valores aproximados para los coeficientes de q(s) en X(s):
den_x = [1.0 21.4258 3.5694 -0.2835 0.0];
% Función de transferencia X(s)
sys_x = tf(num_x, den_x);

% posición cañón: theta
% Valores aproximados para los coeficientes de p(s) acoplado de Theta(s):
num_theta = 10^-5*[2.1195 0.4817];
% Valores aproximados para los coeficientes de p(s) acoplado de Theta(s):
den_theta = [1.0 21.4258 3.5694 -0.2835];
% Función de transferencia Theta(s)
sys_theta = tf(num_theta, den_theta);

% coupling: 
% Valores aproximados para los coeficientes de p(s) en X(s):
num_thetax = 10^-6*2.8902;
% Valores aproximados para los coeficientes de q(s) en X(s):
den_thetax = [1.0 21.4258 3.5694 -0.2835];
% Coupling TF:
sys_thetax = tf(num_thetax, den_thetax); 

% bode: 
% x:
figure;
bode(sys_x); 
title("Bode: control movimiento horizontal carro")
grid on; 
% theta: 
figure; 
bode(sys_theta); 
title("Bode: control de elevación cañón"); 
grid on; 
% coupling:
figure; 
bode(sys_thetax); 
title("Bode: control de acoplamiento señales curzadas"); 
grid on; 
