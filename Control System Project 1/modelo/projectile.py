import numpy as np
import pygame 
import time 
from scipy.integrate import odeint


# ---- Pantalla ----
XMAX = 1000
XMIN = 0
ZMAX = 500
ZMIN = 0
LAMBDA = 0.05       # 20 [pix] = 1 [m]; m/lambda = pix


# ---- Parámetros ----
params = dict()
params["m_x"] = 40000       # masa carro [kg]
params["m_c"] = 4000        # masa cañón [kg]
params["m_b"] = 100         # masa bombarda [kg]
params["b"] = 10**6         # coeficiente de fricción pivote [Nm/(rad/s)]
params["c"] = 10**4         # coeficiente de fricción carro [N/(m/s)]
params["d"] = 1.15          # coeficiente de fricción bombarda [N/(m/s)]
params["g"] = 9.81          # aceleración de gravedad [m/s^2]
params["l_c"] = 6           # largo canón [m]
params["r"] = 0.1           # radio bombarda [m]


# ---- Bomb params ----
QMAX = 6*(10**5)     # carga máxima explosiva [N]
QMIN = 0             # carga mínima explosiva [N]
x = 0                # pos. inicial eje horizontal bombarda [m]
z = 0                # pos. inicial eje vertical bombarda [m]
x_dot = 0            # vel. inicial eje horizontal bombarda [m/s]
z_dot = 0            # vel. inicial eje vertical bombarda [m/s]
r_bomb = 5           # water bomb radius [pix]


# ---- Sampling Period ----
Ts = 0.1


# ---- Clase Bombarda ----

class Projectile(): 

    def __init__(self, wb_init, theta_ref, v_x0, v_z0, r_bomb=r_bomb):
        self.x_state = np.array([0, 0, v_x0, v_z0])  # x_satate = (x, z, x_dot, z_dot)
        self.theta = theta_ref
        self.r_bomb = r_bomb
        self.wb_init = wb_init

    def update_state_projectile(self): 

        # Kinematics of the water bomb:
        def kinematics_model(x, t, params): 

            x_b = x[0]
            z_b = x[1]
            x_b_dot = x[2]
            z_b_dot = x[3]
            m_b = params["m_b"]
            d = params["d"]
            g = (-1)*params["g"]
            x_b_ddot = (1/m_b)*((-1)*d*x_b_dot)
            z_b_ddot = (1/m_b)*((-1)*d*z_b_dot-m_b*g)
            return np.array([x_b_dot, z_b_dot, x_b_ddot, z_b_ddot])

        x_state_0 = self.x_state
        t = np.linspace(0, Ts, 2)
        x_state_b = odeint(kinematics_model, x_state_0, t, 
                        args=(params,))
        
        x_b = x_state_b[-1,0]
        z_b = x_state_b[-1,1]
        x_b_dot = x_state_b[-1,2]
        z_b_dot = x_state_b[-1,3]

        self.x_state = np.array([x_b, z_b, x_b_dot, z_b_dot])

    def check_state_bounds(self, x_pix, z_pix):
        """
        Check bounds in pixel domain: 
        """
        if x_pix >= XMAX or x_pix <= XMIN: 
            self.r_bomb = 0
        if z_pix >= ZMAX or z_pix <= ZMIN: 
            self.r_bomb = 0

    def projectile_pos_display(self, x, z): 
        x_pix = int(x*LAMBDA)
        z_pix = int(z*LAMBDA)
        return x_pix, z_pix
    
    def update_projectile(self): 
        """
        Updates the display of the motion of the water bomb in pixel domain.
        """
        self.update_state_projectile()
        x = self.x_state[0]
        z = self.x_state[1]
        x_b, z_b = self.projectile_pos_display(x, z)
        x_b += self.wb_init[0]
        z_b += self.wb_init[1]
        self.check_state_bounds(x_b, z_b)
        info_projectile = {"[pix]": [x_b, z_b, self.r_bomb]}
        return info_projectile

