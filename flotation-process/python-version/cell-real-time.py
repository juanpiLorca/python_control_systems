import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import json 
import time
import csv 
import keyboard


class Cell(): 
    """
    Flotation process real-time simulator: 
    Real time user inputs given by the keyboard keys. 
        - Use "w" to increase and "s" to decrease manipulated variable valve aperture.
        - Use "d" to increase and "a" to decrease manipulates variavle air flux.
    Run plot-real-time-data script to visualize the behaviour of the states.
    args: 
        >>> idx: index for state-space vectors
        >>> num_classes: number of bubble classes
        >>> sampling_period: sampling period to compute ODEs solutions
    """

    def __init__(self, cell_id: int, i: int, num_classes: int,
                 sampling_period: float, u1Var: float, u2Var: float):
        self.cell_id = cell_id
        self.i = i
        self.K = num_classes
        # Model parameters: 
        with open("params.json") as fjson: 
            params = json.load(fjson)
        self.params = params
        epsis = np.array(params["epsis"])

        # Control: 
        self.Ts = sampling_period
        self.t = 0
        self.x0 = np.concatenate([np.zeros(i), [0.13], epsis / (1 - np.sum(epsis)), [1e-6]])
        self.x = self.x0.copy()
        self.u = np.array([0.3, 0.0013])
        self.u1Var = u1Var
        self.u2Var = u2Var

        # Model control variables: 
        # 1.- Masses: (TODO: masses are not a one dim. array)
        self.m = np.zeros(i)
        self.mfeed = np.zeros(i)
        self.mtail = np.zeros(i)
        self.mTF = np.zeros(i)
        self.mENT = np.zeros(i)
        # 2.- Air-Recovery: 
        self.alpha = 0
        # 3.- Level pulp in cell: 
        self.hp = 0
        # 4.- Air-flow: 
        self.Qair_k = np.zeros(num_classes)

        # Model concentration variables
        self.Qconc = 0 # ???
        self.Qtail = 0
        self.Qpulp_out = 0

        # Data collection: 
        self.time_lapsed = []
        # 1.- Control variables:
        self.ore_mass_data = []
        self.mfeed_data = []
        self.mtail_data = []
        self.mTF_data = []
        self.mENT_data = []
        self.alpha_data = []
        self.hp_data = []
        self.Qair_k_data = []
        # 2.- Manipulated variables: 
        self.valve_aperture = []
        self.air_flux = []
        # 3.- Concentration data:
        self.Qtail_data = []
        self.Qconc_data = []

    def odefun(self, x, t): 
        # State-Space variables: 
        self.m = x[:self.i]
        # h0 = x[self.idx]
        gh = x[self.i+1:self.i+self.K+1]
        dbf = x[self.i+self.K+1]
        # State-Space derivative: 
        DdbfDt = self.params["C"] * dbf**(1 - self.params["n"])

        # Algebraic loop solutions: 
        self.time_lapsed.append(t)
        s = self.model_equations(t, x)
        self.mfeed = s[0:2]
        self.mtail = s[2:4]
        self.mTF = s[4:6]
        self.mENT = s[6:8]
        self.alpha = s[8]
        self.hp = s[9]
        self.Qair_k = s[10:15]
        self.Qconc = s[15]
        self.Qtail = s[16]
        gh_tot = s[17]
        vg_k_pulp_out = s[18:23]
        not_solve = s[23]

        if not not_solve: 
            Qpulp_out = self.Qtail + self.Qconc
            Dh0Dt = (1 / self.params["Acell"]) * (self.params["Qfeed"] - Qpulp_out)
            DghDt = ((1 + gh_tot) / (self.params["Acell"] * self.hp)) * (self.Qair_k - self.params["Acell"]
                        * vg_k_pulp_out * (gh / (1 + gh_tot)) - (self.params["Qfeed"] - Qpulp_out) * gh)
            DghDt = DghDt.reshape((self.K,))
            DmDt = self.mfeed - self.mtail - self.mTF - self.mENT
        
        else: 
            Dh0Dt = 0
            DghDt = np.zeros(self.K)
            DmDt = np.zeros(self.i)

        return np.concatenate([DmDt, [Dh0Dt], DghDt, [DdbfDt]], axis=0)

    def model_equations(self, x): 
        m = x[:self.i]
        h0 = x[self.i]
        gh = x[self.i+1:self.i+self.K+1]
        dbf = x[self.i+self.K+1]

        gh_tot = np.sum(gh)
        epsi0 = gh / (1 + gh_tot)
        uv = self.u[0]  
        Qair_in = self.u[1]  

        gh_tot = np.sum(gh)
        epsi0 = gh / (1 + gh_tot)
        epsi0_tot = np.sum(epsi0)
        hp = h0 / (1 - epsi0_tot)

        density_pulp = (self.params["phi"] * self.params["density_sol"] 
                        + (1 - self.params["phi"]) * self.params["density_water"])
        viscosity_pulp = (self.params["viscosity_water"] 
                          * np.exp(2.5 * self.params["phi"]/(1 - 0.609 * self.params["phi"])))

        hf = self.params["hT"] - hp
        vg_k_pulp_out = ((self.params["g"] * density_pulp * np.array(self.params["db_k_pulp"])**2) 
                         / (18 * viscosity_pulp * (1 - epsi0)**1.39))
        db_int = np.sum(vg_k_pulp_out * epsi0) / np.sum(vg_k_pulp_out * epsi0 / self.params["db_k_pulp"])
        jg = Qair_in/self.params["Acell"]
        vb = self.params["a"] + self.params["b"] * jg + self.params["c"] * jg**2
        K1 = density_pulp * self.params["g"]/ (3 * self.params["CPB"] * viscosity_pulp)
        Qtail = self.params["kv"] * uv * np.sqrt(h0)
        vg_total_out_pulp = np.sum(vg_k_pulp_out * epsi0)

        psi_k = (np.array(self.params["freq_bubble_size"]) 
                 / np.sum(np.array(self.params["freq_bubble_size"])))
        Qair_k = Qair_in * psi_k

        def algebraic_loop(vars): 
            vg_, alpha = vars
            alpha_ = (vg_ - vb) / vg_
            alpha_sat = (self.params["aSAT"] + (self.params["bSAT"] 
                         / (1 + np.exp(-self.params["cSAT"] * (alpha_ - self.params["dSAT"])))))
            tau_f = hf / vg_
            dbf_out = ((self.params["n"] * self.params["C"] * tau_f 
                        + db_int**self.params["n"])**(1/self.params["n"]))
            lambda_out = self.params["klambda"] / dbf_out**2
            Sig = 1 / (1 + np.exp(-self.params["aSF"] * (alpha - self.params["bSF"])))
            Qconc = ((1 - Sig) * (self.params["Acell"] * vg_**2 * lambda_out) 
                     * (1 - alpha_sat) * alpha_sat / K1) + (Sig * self.params["Acell"] 
                        * vg_**2 * lambda_out / 4 * K1)
            DghDt = (((gh_tot + 1) / (self.params["Acell"] * hp)) 
                     * (Qair_k - self.params["Acell"] * vg_k_pulp_out * (gh / (1 + gh_tot)) 
                        - (self.params["Qfeed"] - Qtail - Qconc) * gh))
            epsi = ((1 - Sig) * vg_ * (1 - alpha_sat) * lambda_out / K1) + (Sig * vg_ * lambda_out / (2 * K1))
            eqn1 = (vg_ - (h0 * np.sum(DghDt) - ((gh_tot + 1) 
                    * (self.params["Qfeed"] - Qtail - Qconc)) / self.params["Acell"]) - vg_total_out_pulp)
            eqn2 = alpha - (Qconc / (epsi * Qair_in))
            return [eqn1, eqn2]
        
        init_guess = [0.003, 0.03]
        solution = fsolve(algebraic_loop, init_guess)

        if not np.any(np.isnan(solution)): 
            not_solve = 0
            vg_, alpha = solution
            alpha_ = (vg_ - vb) / vg_
            alpha_sat = (self.params["aSAT"] + (self.params["bSAT"] 
                         / (1 + np.exp(-self.params["cSAT"] * (alpha_ - self.params["dSAT"])))))
            tau_f = hf / vg_
            dbf_out = (self.params["n"] * self.params["C"] * tau_f 
                       + db_int**self.params["n"])**(1/self.params["n"])
            lambda_out = self.params["klambda"] / dbf_out**2
            Sig = 1 / (1 + np.exp(-self.params["aSF"] * (alpha - self.params["bSF"])))
            Qconc = (((1 - Sig) * (self.params["Acell"] * vg_**2 * lambda_out) 
                      * (1 - alpha_sat) * alpha_sat / K1) 
                      + (Sig * self.params["Acell"] * vg_**2 * lambda_out / 4 * K1))
        else:
            not_solve = 1

        if not not_solve:
            Daxial = (jg**1.5) / (np.sqrt(K1 * (np.sqrt(3) - np.pi / 2)) * self.params["Pe"])
            Sb = 6 * vg_ / db_int
            vset = (self.params["g"] * (self.params["density_sol"] - self.params["density_water"]) 
                    * np.array(self.params["dp"])**2 * (1 - self.params["phi"])**4.65 
                    / (18 * viscosity_pulp * 3))
            Rf = (((1 - Sig) * (((alpha_sat * (1 - alpha_sat) * vg_) / vset)**(self.params["f"] / 2)) 
                   * ((db_int / dbf)**self.params["f"])) + Sig * ((vg_ / vset)**(self.params["f"] / 2)) 
                   * ((db_int / dbf)**self.params["f"]))
            ENT = (((1 - Sig) * np.exp(-(vset**1.5) * hf / (Daxial * np.sqrt(vg_ * (1 - alpha_sat))))) 
                   + Sig * np.exp(-2 * (vset**1.5) * hf / (Daxial * np.sqrt(vg_))))
            k = np.array(self.params["P"]) * Sb
        
        Vk_gas = (epsi0 / (1 - epsi0)) * h0 * self.params["Acell"]
        Vgas = np.sum(Vk_gas)
        Vpulp = h0 * self.params["Acell"] + Vgas
        Ctail = m / Vpulp

        # Masses: Mass Balance 
        mfeed = np.array(self.params["Cfeed"]) * np.array(self.params["Qfeed"])
        mtail = Ctail * Qtail
        mTF = self.params["Vcell"] * k * Rf * Ctail
        mENT = Qconc * ENT * Ctail

        # Data collection: 
        # 1.- 
        self.ore_mass_data.append(m)
        self.mfeed_data.append(mfeed[0])
        self.mtail_data.append(mtail[0])
        self.mTF_data.append(mTF[0])
        self.mENT_data.append(mENT[0])
        self.alpha_data.append(alpha)
        self.hp_data.append(hp)
        self.Qair_k_data.append(np.sum(Qair_k))
        # 2.- 
        self.valve_aperture.append(self.u[0])
        self.air_flux.append(self.u[1])
        # 3.-
        self.Qtail_data.append(Qtail)
        self.Qconc_data.append(Qconc)

        return np.concatenate([mfeed, mtail, mTF, mENT, [alpha], [hp], Qair_k, 
                               [Qconc], [Qtail], [gh_tot], vg_k_pulp_out, [not_solve]], 
                               axis=0)
    
    def actuators(self, e):
        u1 = self.u[0]
        u2 = self.u[1]
        if e.event_type == keyboard.KEY_DOWN:
            if e.name == "w": 
                u1 += self.u1Var
            if e.name == "s": 
                u1 -= self.u1Var
            if e.name == "d": 
                u2 += self.u2Var
            if e.name == "a": 
                u2 -= self.u2Var

        u1, u2 = self.limits(u1, u2)
        self.u = np.array([u1, u2])

    def limits(self, u1, u2): 
        if u1 >= 1.0: 
            u1 = 1.0 
        if u1 <= 0.0: 
            u1 = 0.0
        if u2 >= 0.002: 
            u2 = 0.002 
        if u2 <= 0.001: 
            u2 = 0.001
        return u1, u2

    def simulate(self): 
        self.x0 = self.x 
        t = np.linspace(0, self.Ts, 2)
        # Perform integration using Fortran's LSODA (Adams & BDF methods)
        x = odeint(self.odefun, self.x0, t) 
        m = x[-1,:self.idx]
        h0 = x[-1,self.idx]
        gh = x[-1,self.idx+1:self.idx+self.K+1]
        dbf = x[-1,self.idx+self.K+1]
        self.x = np.concatenate([m, [h0], gh, [dbf]], axis=0)
        self.time_lapsed.append(self.t)
        self.t += self.Ts
        return self.x
    

fieldnames = ["time", "y_value"]
with open("data.csv", "w") as csv_file:
    csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    csv_writer.writeheader() 


def main():

    i = 2
    K = 5
    Ts = 1.0
    valve_aperture = 5.0e-4
    air_flux = 1.0e-6
    cell = Cell(idx=i, 
                num_classes=K, 
                sampling_period=Ts, 
                u1Var=valve_aperture, 
                u2Var=air_flux)

    while True: 
        keyboard.hook(cell.actuators)
        x = cell.simulate()
        y = x[0]
        t = cell.t
        with open("data.csv", "a") as csv_file: 
            csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            info = {
                "time": t, 
                "y_value": y
            }      
            csv_writer.writerow(info)
            print(f"Opening valve: {round(cell.u[0], 5)}, Air flux: {round(cell.u[1], 5)}")

        time.sleep(0.5)

if __name__ == "__main__": 
    main()
    
