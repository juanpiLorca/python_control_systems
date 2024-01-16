% Not using Miguel Torres coefficients in the TF: 

% Sampling period: 
Ts = 0.1;

% x: 
num_x = [5.76*10^5 1.2*10^5 0.0]; 
den_x = [2.534*10^10 5.338*10^11 1.2*10^11 0.0 0.0]; 
sysx = tf(num_x, den_x); 
sysDx = c2d(sysx, Ts);

% theta: 
num_theta = [5.28*10^5 1.2*10^7]; 
den_theta = [2.534*10^10 5.338*10^11 1.2*10^11 0.0]; 
systheta = tf(num_theta, den_theta); 
sysDtheta = c2d(systheta, Ts); 

% Bode: x-position 
figure; 
bode(sysDx); 
title("bode: control de posición horizontal carro"); 
grid on; 

% Bode: theta-elevation; 
figure; 
bode(sysDtheta); 
title("bode: control de elevación cañón"); 
grid on; 