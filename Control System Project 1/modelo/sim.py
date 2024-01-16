from model import Model
from fire import Fire
from projectile import Projectile
import numpy as np
import pygame 
import time

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


# ---- Pantalla ---- [pix]
XMAX = 1000
ZMAX = 500
LAMBDA = 0.05       # 20 [pix] = 1 [m]; m/lambda = pixs
screen = None
tw0 = 0
click_radius = 4    


# ---- Carro y cañón dim ---- [pix]
CART_WIDTH = 30
CART_HEIGHT = 20
WHEELS_DIAMETER = 10


# ---- Variables Cañón-Carro ----
x = 0 
theta = 0
x_dot = 0
theta_dot = 0


# ---- Control Variables ----
mv_cart = 0
mv_cannon = (-1)*params["m_c"]*params["g"]*params["l_c"]/2*np.cos(theta)
mv_load = 0


class ModelDisplayPixel():

    def __init__(self): 
        self.screen = None
        self.model = Model(x, theta, x_dot, theta_dot, 
                           mv_cart, mv_cannon, mv_load)
        self.cart_cm = None

        # Input clicks register: 
        self.input_clicks = []
        
        # Fire register: 
        self.fire_reg = []
        self.handle_fire = False

        # Projectile register: water bomb register
        self.init_pos_wb = None
        self.wb_reg = []
        
    def init_display(self): 
        pygame.init()
        self.screen = pygame.display.set_mode((XMAX, ZMAX))
        pygame.display.set_caption("cannon 23 sim")
        pygame.key.set_repeat(1,50)     # debouncer
        self.screen.fill((0,0,0))
        pygame.display.flip()

    def x_display_pix(self, x, z=0):    # cart's center of mass
        x_pix = int(XMAX/2 + x*LAMBDA)
        z_pix = int(ZMAX/2 + z*LAMBDA) - WHEELS_DIAMETER + 2
        x_pix_cm = int(XMAX/2 + CART_WIDTH/2 + x*LAMBDA)
        z_pix_cm = int(ZMAX/2 + CART_HEIGHT/2 + z*LAMBDA) - WHEELS_DIAMETER + 2
        return x_pix, z_pix, x_pix_cm, z_pix_cm
    
    def cannon_display(self, x, z, theta, l_c=params["l_c"]):
        x_pix = x + l_c*np.cos(theta)*5
        z_pix = z + l_c*np.sin(theta)*5
        return x_pix, z_pix

    def model_display(self): 
        info = self.model.update_model()
        x = info["x"]
        theta = info["theta"]
        load = info["load"] 
        x_ref = info["x_ref"]
        theta_ref = info["theta_ref"]
        x_cart, z_cart, x_cm, z_cm = self.x_display_pix(x)
        self.cart_cm = (x_cm, z_cm)
        x_cann, z_cann = self.cannon_display(x_cm, z_cm, theta)
        self.init_pos_wb = (x_cann, z_cann)

        self.screen.fill((0,0,0))
        pygame.draw.rect(self.screen, (100,40,0), 
                         pygame.Rect(0, ZMAX/2 + CART_HEIGHT, XMAX, ZMAX/2))
        # Cart display:
        pygame.draw.rect(self.screen, (130,130,130), 
                         pygame.Rect(x_cart, z_cart, CART_WIDTH, CART_HEIGHT))
        # Cart's wheels: left
        x_wl = x_cart + WHEELS_DIAMETER/2
        z_wl = ZMAX/2 + (CART_HEIGHT - WHEELS_DIAMETER/2) 
        pygame.draw.circle(self.screen, (100, 100, 100), 
                           (x_wl, z_wl), WHEELS_DIAMETER/2)
        # Cart's wheels: right
        x_wr = x_cart + (CART_WIDTH - WHEELS_DIAMETER/2)
        z_wr = ZMAX/2 + (CART_HEIGHT - WHEELS_DIAMETER/2) 
        pygame.draw.circle(self.screen, (100, 100, 100), 
                           (x_wr, z_wr), WHEELS_DIAMETER/2)
        
        # Cannon display: 
        pygame.draw.line(self.screen, (100, 100, 100), start_pos=(x_cm, z_cm), 
                     end_pos=(x_cann, z_cann), width=int(4))
        
        if self.handle_fire:
            if self.model.auto_shoot == 1: 
                self.gen_auto_projectile()
                print("Auto shot done!")

        cart_str = f"Cart x-position: {round(x)} [m]" 
        ref_cart_str = f"Cart x-position reference: {round(x_ref)} [m]"
        angle = np.degrees(theta)
        cannon_str = f"Cannon theta-angle: {round(angle)}° [degrees]"
        angle_ref = np.degrees(theta_ref)
        ref_cannon_str = f"Cannon theta-angle reference: {round(angle_ref)}° [degrees]"
        load_str = f"Load: {load} [N]"
        font = pygame.font.SysFont("Arial", 20)
        self.screen.blit(font.render(cart_str, True, (255,255,255)), (25,25))
        self.screen.blit(font.render(ref_cart_str, True, (255,255,255)), (25,50))
        self.screen.blit(font.render(cannon_str, True, (255,255,255)), (25,75))
        self.screen.blit(font.render(ref_cannon_str, True, (255,255,255)), (25,100))
        self.screen.blit(font.render(load_str, True, (255,255,255)), (25,125))
    
    def gen_auto_projectile(self):
        theta_ref = self.model.shoot_projectile()["theta_init"]
        v_x0 = self.model.shoot_projectile()["v_x0"]
        v_z0 = self.model.shoot_projectile()["v_z0"]
        projectile = Projectile(self.init_pos_wb, theta_ref, v_x0, v_z0)
        self.wb_reg.append(projectile)

    def display_click_input(self): 
        if len(self.input_clicks) > 0: 
            for click in self.input_clicks:
                x_pix = click[0]
                z_pix = click[1]
                pygame.draw.circle(self.screen, (255,255,255,255), 
                                    (x_pix, z_pix), click_radius)

    def fire_display(self): 
        if len(self.fire_reg) > 0: 
            for fire in self.fire_reg: 
                info_pix = fire.fire_pos["[pix]"]
                x_pix, z_pix, r_pix = info_pix[0], info_pix[1], info_pix[2]
                pygame.draw.circle(self.screen, (255,100,10), 
                                    (x_pix, z_pix), r_pix)
                pygame.draw.circle(self.screen, (255,255,0), 
                                    (x_pix, z_pix), r_pix*0.7)
                
    def squared_range(self, x_b, z_b): 
        """
        Checks if the projectile has reached the fire.
        """
        hit = False
        for fire in self.fire_reg:
            x, z, r = fire.fire_pos["[pix]"][0], fire.fire_pos["[pix]"][1], fire.fire_pos["[pix]"][2]
            x_range = (x - r, x + r)
            z_range = (z - r, z + r)
            if (x_range[0] <= x_b <= x_range[1]) and (z_range[0] <= z_b <= z_range[1]): 
                hit = True
                self.handle_fire = False
                self.model.auto_shoot = 0
            return hit
                
    def projectile_display(self): 
        if len(self.wb_reg) > 0: 
            for water_bomb in self.wb_reg: 
                info_pix = water_bomb.update_projectile()
                info_pix = info_pix["[pix]"]
                x_pix, z_pix, r_pix = info_pix[0], info_pix[1], info_pix[2]
                pygame.draw.circle(self.screen, (0,0,255), (x_pix, z_pix), r_pix)
                if self.squared_range(x_pix, z_pix): 
                    self.fire_reg.pop()
    
    def handle_keyboard(self): 
        for event in pygame.event.get(): 
            if event.type == pygame.QUIT: 
                    pygame.quit();

            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    pygame.quit();
                
                # Model:
                if event.key == pygame.K_RIGHT or event.key == pygame.K_LEFT:
                    self.model.keyboard_logic(event.key)
                if event.key == pygame.K_UP or event.key == pygame.K_DOWN:
                    self.model.keyboard_logic(event.key)
                if event.key == pygame.K_1 or event.key == pygame.K_2: 
                    self.model.keyboard_logic(event.key)
                if event.key == pygame.K_s: 
                    self.model.keyboard_logic(event.key)
                if event.key == pygame.K_a: 
                    self.model.keyboard_logic(event.key)  

                # Fire:
                if event.key == pygame.K_f:
                    fire = Fire()
                    fire.random_fire()
                    if len(self.input_clicks) >= 1: 
                        self.input_clicks.pop()
                    if len(self.fire_reg) < 1:
                        self.fire_reg.append(fire)
                        self.handle_fire = True
                        x_ref = fire.fire_pos["[pix]"][0] # [pix]
                        z_ref = fire.fire_pos["[pix]"][1] # [pix]
                        # x_reference for the cart-cannon model will be at 2 [km] from the fire
                        if x_ref >= self.cart_cm[0]: 
                            x_ref -= 2000*LAMBDA
                            x_angle = 2000*LAMBDA - CART_WIDTH/2
                        if x_ref <= self.cart_cm[0]: 
                            x_ref += 2000*LAMBDA
                            x_angle = (-1)*2000*LAMBDA - CART_WIDTH/2
                        x_reference = (x_ref - XMAX/2)/LAMBDA 
                        self.model.x_reference = x_reference
                        # theta_reference: 
                        z_angle = z_ref - self.cart_cm[1]
                        theta_reference = np.arctan2(z_angle, x_angle)
                        self.model.theta_reference = theta_reference

                # Projectiles: water bombs
                if event.key == pygame.K_SPACE:
                    if self.model.mv_load > 0: 
                        self.model.keyboard_logic(event.key)
                        theta_ref = self.model.shoot_projectile()["theta_init"]
                        v_x0 = self.model.shoot_projectile()["v_x0"]
                        v_z0 = self.model.shoot_projectile()["v_z0"]
                        projectile = Projectile(self.init_pos_wb, theta_ref, v_x0, v_z0)
                        self.wb_reg.append(projectile)
                    else: 
                        print("There's no load to fire water bombs!") 
                        print("Please try again!")
            
            # Graphic point: 
            if event.type == pygame.MOUSEBUTTONDOWN: 
                if event.button == 1: 
                    pos = pygame.mouse.get_pos()
                    ref_pix = (pos[0], pos[1])
                    ref_meters = (round(pos[0]/LAMBDA), round(pos[1]/LAMBDA))
                    print(f"You've clicked on the position {ref_pix} [pix]; {ref_meters} [m]")
                    if not self.handle_fire: 
                        x_ref = ref_pix[0]
                        z_ref = ref_pix[1]
                        # x_reference for the cart-cannon model will be at 2 [km] from the point graphed
                        if x_ref >= self.cart_cm[0]: 
                            x_ref -= 2000*LAMBDA
                            x_angle = 2000*LAMBDA - CART_WIDTH/2
                        if x_ref <= self.cart_cm[0]: 
                            x_ref += 2000*LAMBDA
                            x_angle = (-1)*2000*LAMBDA - CART_WIDTH/2
                        x_reference = (x_ref - XMAX/2)/LAMBDA
                        self.model.x_reference = x_reference
                        # theta_reference: 
                        z_angle = z_ref - self.cart_cm[1]
                        theta_reference = np.arctan2(z_angle, x_angle)
                        self.model.theta_reference = theta_reference
                        if len(self.input_clicks) >= 1: 
                            self.input_clicks.pop()
                        self.input_clicks.append(ref_pix)
                    else: 
                        print("You need to put down the fire first!")
                        print("Try again!")

    def display_scenario(self): 
        self.model_display()
        self.display_click_input()
        self.fire_display()
        self.projectile_display()
        pygame.display.flip()                    


def main():
    tw0 = time.time()
    model_display = ModelDisplayPixel()
    model_display.init_display()
    while True: 
        model_display.handle_keyboard()
        model_display.display_scenario()
        model_display.model.store_data(tw0)

if __name__ == "__main__": 
    main()
    