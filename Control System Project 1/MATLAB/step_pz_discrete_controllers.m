% Discrete controllers for x-position adn theta-elevation: 

% Sampling rate: 
Ts = 0.1;
z = tf('z',Ts);

% x-position: 
sysX = (87226.3*z^2-155000*z+85250)/(Ts*z*(z-1));
figure; 
pzplot(sysX); 
title("Polos y ceros controlador posición horizontal carro");
grid on;

% theta-elevation: 3613.5 15840 330
sysTheta = (52470*z^2-82500+49500)/(Ts*z*(z-1));
figure;
pzplot(sysTheta); 
title("Polos y ceros controlador ángulo de elevación cañón");
grid on;
