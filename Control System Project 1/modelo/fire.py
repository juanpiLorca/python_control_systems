import numpy as np
import pygame 
import time 


# ---- Pantalla ----
XMAX = 1000
XMIN = 0
ZMAX = 500
ZMIN = 0
LAMBDA = 0.05       # 20 [pix] = 1 [m]; m/lambda = pix


# ---- Fuego ----
RADIUS = 5


# ---- Clase Generadora de Fuego ----

class Fire(): 

    def __init__(self):
        self.min_radius_pix = RADIUS
        self.min_radius_m = RADIUS/LAMBDA
        self.fire_pos = {}

    def random_fire(self):
        ran_fire = (np.array([10000, 2500, 2000])*2*
                    (np.random.rand(3)+np.array([-0.5, -0.5, -0.5])))
        x_pix = round(ran_fire[0]*LAMBDA + XMAX/2)
        if x_pix >= XMAX: 
            x_pix = XMAX
        if x_pix <= XMIN: 
            x_pix = XMIN
        self.x = x_pix
        z_pix = round(ran_fire[1]*LAMBDA + ZMAX/2)
        if z_pix >= ZMAX: 
            z_pix = ZMAX
        if z_pix <= ZMIN: 
            z_pix = ZMIN
        r_pix = round(ran_fire[2]*LAMBDA)
        if r_pix < self.min_radius_pix: 
            r_pix = self.min_radius_pix  
        self.fire_pos["[pix]"] = [x_pix, z_pix, r_pix]
        
    