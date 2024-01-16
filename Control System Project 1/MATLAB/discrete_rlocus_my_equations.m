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

% poles & zeros: 
% x-position:
figure; 
pzplot(sysDx); 
title("Polos y ceros sistema de control de posición horizontal")
grid on; 
% theta-elevation: 
figure; 
pzplot(sysDtheta); 
title("Polos y ceros sistema de control de elevación cañón")
grid on; 

% rlocus: x-position 
figure; 
rlocus(sysDx); 
title("rlocus discrete: control de posición horizontal carro"); 
grid on; 

% rlocus: theta-elevation; 
figure; 
rlocus(sysDtheta); 
title("rlocus discrete: control de elevación cañón"); 
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
sysDTheta = feedback(sysdtheta, 1); 
figure; 
step(sysDTheta); 
title("Step response: closed loop theta-elevation")
grid on; 
