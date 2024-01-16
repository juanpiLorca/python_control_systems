import numpy as np
import pygame 
import time
import os
from scipy.integrate import odeint


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


# ---- Pantalla ----
XMAX = 1000
ZMAX = 500
LAMBDA = 0.05       # 20 [pix] = 1 [m]; m/lambda = pix


# ---- Carro y cañón dim ----
CART_WIDTH = 30
CART_HEIGHT = 20


# ---- Water bomb params ----
QMAX = 6*(10**5)     # carga máxima explosiva [N]
QMIN = 0             # carga mínima explosiva [N]


# ---- Boundaries ----
X_MAX_CART = XMAX/2             # [pixel]
X_MIN_CART = (-1)*(XMAX/2)      # [pixel]

F_XMAX = 30000          # fuerza máxima motor [N]
F_XMIN = -30000         # fuerza mínima motor [N] (para hacer que retroceda)
T_THETAMAX = 235440     # torque máximo motor pivote [N/m]
T_THETAMIN = -235440    # troque mínimo motor pivote [N/m] (para lograr un grio anti-horario)


# ---- Sampling Period ----
Ts = 0.1


# ---- Logger variables ----
k_data = 0                  # Índice de la fila donde se almacenará la data
k_data_samples = 5000       # Número de muestras a almacenar. 5000*0.02 = 100 [s] sim.
record_size = 7             # Cantidad de datos a almacenar en una fila: time+state_vector+referece+manipulated_variable = 1+2+2+2
sampling_interval = 0.02    # [s]
sim_time = 0                # [s]


# ---- Controller params: PID ----
# x-system:
K_x = 1.55e5
Kp_x = 1.0*K_x        # proportional gain
Ki_x = 0.275*K_x      # integral proportional gain
Kd_x = 0.55*K_x       # derivative proportional gain
# theta-system:
K_theta = 1.7e5
Kp_theta = 1.0*K_theta      # proportional gain
Ki_theta = 0.8*K_theta      # integral proportional gain
Kd_theta = 0.3*K_theta      # derivative proportional gain


# ---- Noise params ----
A = 10e-3
ST_DV = 1


class Model(): 

    def __init__(self, x, theta, x_dot, theta_dot, 
                 mv_cart, mv_cannon, mv_load):
        # state model atributes:
        self.x_state = np.array([x, theta, x_dot, theta_dot])
        self.mv_cart = mv_cart
        self.mv_cannon = mv_cannon
        self.mv_load = mv_load
        # refrences atributes:
        self.x_reference = self.x_state[0]
        self.theta_reference = self.x_state[1]
        # error atributes:
        self.x_error = 0
        self.x_error_old = 0
        self.x_error_old2 = 0
        self.theta_error = 0
        self.theta_error_old = 0
        self.theta_error_old2 = 0
        # automatic mode atributes:
        self.auto = 0
        self.auto_shoot = 0
        # storing data atributes:
        self.data = [[],[],[],[],[],[],[]]  # len(self.data) = record_size = 7
        self.k_data = k_data
        self.sim_time = sim_time

    def update_state_model(self): 

        # Kinematics of the cart-cannon model: 
        def kinematics_model(x_state, t, u, params): 
            x = x_state[0]
            theta = x_state[1]
            x_dot = x_state[2]
            theta_dot = x_state[3]
            F_x = u[0]
            F_c = u[1]
            m_x = params["m_x"]
            m_c = params["m_c"]
            m_T = m_x + m_c
            l_c = params["l_c"]
            g = (-1)*params["g"]             
            b = params["b"]
            c = params["c"]

            x_ddot = ((4*l_c*(F_x-c*x_dot+(l_c/2)*m_c*theta_dot**2*np.cos(theta)) + 
                        6*np.sin(theta)*(F_c-b*theta_dot-g*(l_c/2)*m_c*np.cos(theta)))/
                        (l_c*(3*m_c*np.cos(theta)**2+m_c+4*m_x)))
            theta_ddot = ((6*(l_c*m_c*np.sin(theta)*(F_x-c*x_dot+(l_c/2)*m_c*theta_dot**2*np.cos(theta)) + 
                            4*m_T*(F_c-b*theta_dot-g*(l_c/2)*m_c*np.cos(theta))))/
                            (l_c**2*m_c*(3*m_c*np.cos(theta)**2+m_c+4*m_x)))
            return np.array([x_dot, theta_dot, x_ddot, theta_ddot])
        
        x_state_0 = self.x_state
        mv_u = [self.mv_cart, self.mv_cannon]
        t = np.linspace(0, Ts, 2)
        x_state = odeint(kinematics_model, x_state_0, t, args=(mv_u, params))
        x = x_state[-1,0]
        x = self.check_state_bounds(x)
        theta = x_state[-1,1]
        x_dot = x_state[-1,2]
        theta_dot = x_state[-1,3]
        self.x_state = np.array([x, theta, x_dot, theta_dot])

    def check_state_bounds(self, x):  
        """
        Returns boundaries for x in pixel domain.
        """
        x_pix = int(x*LAMBDA)
        if x_pix >= X_MAX_CART: 
            x = X_MAX_CART/LAMBDA
        if x_pix <= X_MIN_CART: 
            x = X_MIN_CART/LAMBDA
        return x

    def update_model(self): 
        self.update_state_model()

        if self.auto == 1: 
            self.x_controller()
            self.theta_controller()

        x = self.x_state[0]
        theta = self.x_state[1]

        info_model = {"x": x, "theta": theta, "load": self.mv_load, 
                      "x_ref": self.x_reference, "theta_ref": self.theta_reference}
        return info_model

    def e_stopping(self): 
        x = self.x_state[0] 
        theta = self.x_state[1]
        self.mv_cart = 0    
        self.mv_cannon = 0
        self.x_state = np.array([x, theta, 0, 0]) 

    def x_controller(self): 
        self.x_error = self.x_reference - self.x_state[0]
        self.mv_cart += (Kp_x*(self.x_error-self.x_error_old) + Ki_x*Ts*self.x_error 
                         + Kd_x*(self.x_error-2*self.x_error_old+self.x_error_old2)/Ts)
        self.x_error_old2 = self.x_error_old
        self.x_error_old = self.x_error
        # Saturation curve: 
        if self.mv_cart >= F_XMAX:
            self.mv_cart = F_XMAX
        if self.mv_cart <= F_XMIN: 
            self.mv_cart = F_XMIN

    def theta_controller(self): 
        self.theta_error = self.theta_reference - self.x_state[1] # + A*np.random.normal(0,ST_DV)
        self.mv_cannon += (Kp_theta*(self.theta_error-self.theta_error_old) + Ki_theta*Ts*self.theta_error 
                        + Kd_theta*(self.theta_error-2*self.theta_error_old+self.theta_error_old2)/Ts)
        self.theta_error_old2 = self.theta_error_old
        self.theta_error_old = self.theta_error
        # Saturation curve: 
        if self.mv_cannon >= T_THETAMAX:
            self.mv_cannon = T_THETAMAX
        if self.mv_cannon <= T_THETAMIN: 
            self.mv_cannon = T_THETAMIN
        # Reaching reference auto shooting:
        if 0 < self.x_error < 0.01:
            if 0 < self.theta_error < 0.1: 
                self.auto_shoot = 1
            else: 
                self.auto_shoot = 0

    def shoot_projectile(self): 
        theta = self.x_state[1]
        Q_out = self.mv_load
        m_b = params["m_b"]
        v_x0 = (Q_out*Ts*np.cos(theta))/m_b
        v_z0 = (Q_out*Ts*np.sin(theta))/m_b
        init_wb = {"theta_init": theta, "v_x0": v_x0, "v_z0": v_z0}
        return init_wb

    def keyboard_logic(self, key): 
        if key == pygame.K_RIGHT:
            if self.auto == 0:
                if self.mv_cart == F_XMAX: 
                    self.mv_cart = F_XMAX
                    print(f"Engine operating at it's maximum {self.mv_cart} [N]!")
                else: 
                    self.mv_cart += 500
                    print(f"Engine operating at {self.mv_cart} [N]") 
            else: 
                self.x_reference += 0.05*(20*10**3) # 5% of the maximun distance
                print(f"x-reference set at: {self.x_reference} [m]")

        if key == pygame.K_LEFT:
            if self.auto == 0:
                if self.mv_cart == F_XMIN: 
                    self.mv_cart = F_XMIN
                    print(f"Engine operating at it's mininum {self.mv_cart} [N]!")
                else: 
                    self.mv_cart -= 500
                    print(f"Engine operating at {self.mv_cart} [N]")
            else: 
                self.x_reference -= 0.05*(20*10**3) # 5% of the maximun distance
                print(f"x-reference set at: {self.x_reference} [m]")

        if key == pygame.K_DOWN:
            if self.auto == 0: 
                if self.mv_cannon == T_THETAMAX:
                    self.mv_cannon = T_THETAMAX
                    print(f"Pivot torque operating at it's maximum {self.mv_cart} [N/m]!")
                else: 
                    self.mv_cannon += 10**3
                    print(f"Pivote torque operating at {self.mv_cannon} [N/m]")
            else: 
                self.theta_reference += 0.05*(2*np.pi)
                print(f"theta-reference set at: {np.degrees(self.theta_reference)}° [degrees]")

        if key == pygame.K_UP:
            if self.auto == 0:
                if self.mv_cannon == T_THETAMIN: 
                    self.mv_cannon = T_THETAMIN
                    print(f"Pivot torque operating at it's mininum {self.mv_cannon} [N/m]!")
                else: 
                    self.mv_cannon -= 10**3
                    print(f"Pivote torque operating at {self.mv_cannon} [N/m]")
            else: 
                self.theta_reference -= 0.05*(2*np.pi)
                print(f"theta-reference set at: {np.degrees(self.theta_reference)}° [degrees]")

        if key == pygame.K_1:
            if self.auto == 0:
                if self.mv_load == QMAX: 
                    self.mv_load = QMAX
                    print(f"Maximum load reach: {self.mv_load} [N]")
                else: 
                    self.mv_load += 2*10**4
                    print(f"Explosive load: {self.mv_load} [N]")

        if key == pygame.K_2: 
            if self.auto == 0:
                if self.mv_load == QMIN: 
                    self.mv_load = QMIN
                    print(f"Maninum load reach: {self.mv_load} [N]")
                else: 
                    self.mv_load -= 2*10**4
                    print(f"Explosive load: {self.mv_load} [N]")

        if key == pygame.K_SPACE: 
            if self.auto == 0:
                self.shoot_projectile()

        if key == pygame.K_s:
            print("¡e-stop!")
            self.e_stopping()
            self.auto = 0

        if key == pygame.K_a: 
            if self.auto == 0:
                self.auto = 1
                print("\nControl mode is set to [automatic]")
                print(f"You're at {round(self.x_state[0])} [m], {round(self.x_state[1]*180/np.pi)}° [degrees]")
                print(f"Reference is {round(self.x_reference)} [m], {round(self.theta_reference*180/np.pi)}° [degrees]")
                # We set the mv_load (Q) for further calculations:
                self.mv_cart, self.mv_load = 0, 5*10**5
            else: 
                self.auto = 0
                print("\nControl mode is set to [manual]")
                print("\nCart's engine power shuted down: 0 [N]")
                self.mv_cart = 0
                cannon_torque = self.mv_cannon
                self.mv_cannon = cannon_torque

    def store_data(self, tw0):
        """
        Stores automatic simulation data in a .npy file. 
        The sampling period for each sample stored in will be T = 0.02 [s]. 
        Attempts to save up 100 seconds of simulation.
        """
        current_time = time.time() - tw0
        if current_time - self.sim_time >= sampling_interval:
            self.data[0].append(time.time() - tw0)          # storing time: t
            self.data[1].append(self.x_state[0])            # storing horizontal position: x
            self.data[2].append(self.x_state[1])            # storing cannon elevation: theta
            self.data[3].append(self.x_reference)           # storing horizontal pos. ref: x_ref
            self.data[4].append(self.theta_reference)       # storing cannon elevation ref: theta_ref
            self.data[5].append(self.mv_cart)               # storing mv_cart: mv_x
            self.data[6].append(self.mv_cannon)             # storing mv_cannon: mv_c
            self.k_data += 1
            self.sim_time = current_time
        if self.k_data == k_data_samples: 
            data = np.array(self.data)
            np.save("simdata", data)
            print("Data saved into simdata.npy file")
            self.k_data = 0
