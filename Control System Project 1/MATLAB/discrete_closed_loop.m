% Gráficos de acuerdo al modelo corregido por el profesor: Miguel Torres
Ts = 0.1;
% posición carro: x
% Valores aproximados para los coeficientes de p(s) en X(s):
num_x = 10^-5*[2.3121 48.1696 -2.8353];
% Valores aproximados para los coeficientes de q(s) en X(s):
den_x = [1.0 21.4258 3.5694 -0.2835 0.0];
% Función de transferencia X(s)
sys_x = tf(num_x, den_x);
sysDx = c2d(sys_x, Ts);

% posición cañón: theta
% Valores aproximados para los coeficientes de p(s) en X(s):
num_thetax = 10^-6*2.8902;
% Valores aproximados para los coeficientes de q(s) en X(s):
den_thetax = [1.0 21.4258 3.5694 -0.2835];
% Función de transferencia Theta(s)
sys_theta = tf(num_theta, den_theta);
sysDtheta = c2d(sys_theta, Ts); 

% Valores aproximados para los coeficientes de p(s) coupling:
num_theta = 10^-5*[2.1195 0.4817];
% Valores aproximados para los coeficientes de p(s) coupling:
den_theta = [1.0 21.4258 3.5694 -0.2835];
% coupling TF:
sys_thetax = tf(num_thetax, den_thetax); 
sysDxtheta = c2d(sys_thetax, Ts);


% rlocus: x-position 
figure; 
rlocus(sysDx); 
title("rlocus: discrete x-position"); 
grid on; 

% rlocus: theta-elevation; 
figure; 
rlocus(sysDtheta); 
title("rlocus: discrete theta-elevation"); 
grid on; 

% we'll try different critical gains: 
K_x = 3e+01;
K_theta = 1e+01; 

% x-position
sysdx = series(sysDx, K_x); 
sysDX = feedback(sysdx, 1); 
figure; 
step(sysDX); 
title("Step response: closed loop x-position")
grid on; 

% theta-elevation:
sysdtheta = series(sysDtheta, K_theta);
sysDTheta = feedback(systhetad, 1); 
figure; 
step(sysDTheta); 
title("Step response: closed loop theta-elevation")
grid on; 
