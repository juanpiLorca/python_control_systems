import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import json 
import time

class Cell(): 
    """
    Flotation process simulator: 
    Pre-defined user inputs to stimulate the plant in different states. 
    args: 
        >>> cell_id: id of a cell in a line
        >>> i: mineral type.
        >>> num_classes: bubbles size classes.
        >>> time_span: time span for the simulation. 
        >>> u_time_span: manipulated variables for the simulation. 
        >>> u_time_interval: time interval for each change in mv.
    """

    def __init__(self, cell_id: int,i: int, num_classes: int, 
                 time_span: list,u_time_span: list, u_time_interval: int):
        self.id = cell_id
        self.i = i
        self.K = num_classes

        # Model parameters: 
        with open("python-version/params.json") as fjson: 
            params = json.load(fjson)
        self.params = params

        self.tspan = time_span
        self.y0 = np.concatenate((np.zeros(i), [0.13], params["epsis"] / 
                                  (1 - np.sum(params["epsis"])), [1e-6]))
        self.u_tspan = u_time_span
        self.u = np.array(u_time_span[0])

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

        # For user code-interaction
        self.idx_mv = 0
        self.mv_interval = u_time_interval
        self.iter_counter = 0 
        
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

    def update_model(self):

        def odefun(t, x): 
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
        
        ivp_sol = solve_ivp(odefun, t_span=self.tspan, y0=self.y0, method="RK45")
        return ivp_sol

    def model_equations(self, t, x): 
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

        # Maniupalted Variable modifications
        # Works only for t_span = (0,4e4) and 4 inputs
        self.iter_counter += 1
        if t > (1 + self.idx_mv) * self.mv_interval: 
            self.idx_mv += 1
            print(self.idx_mv)
            print(self.u_tspan)
            self.u = np.array(self.u_tspan[self.idx_mv])

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
    


def plot_control_variables(cell: Cell): 
    fig, axs = plt.subplots(nrows=3, ncols=3)
    axs[0,0].plot(cell.ore_mass_data[0], label="Ore") 
    axs[0,0].plot(cell.ore_mass_data[1], label="No ore") 
    axs[0,0].set_title("Masa total de cada mineral")
    axs[0,0].legend()

    axs[0,1].plot(cell.mfeed_data, label="Mfeed")
    axs[0,1].plot(cell.mtail_data, label="Mtail")
    axs[0,1].set_title("Mfeed y Mtail")
    axs[0,1].legend()

    axs[0,2].plot(cell.mTF_data, label="Mtf")
    axs[0,2].plot(cell.mENT_data, label="Ment")
    axs[0,2].set_title('Mtf y Ment')
    axs[0,2].legend()

    axs[1,0].plot(cell.hp_data, label="hp")
    axs[1,0].set_title('Altura de la pulpa en la celda')
    axs[1,0].legend()

    axs[1,1].plot(cell.alpha_data, label="alpha")
    axs[1,1].set_title('Recuperacion de aire')
    axs[1,1].legend()

    axs[1,2].plot(cell.Qtail_data, label="Qtail")
    axs[1,2].set_title('Flujo de cola')
    axs[1,2].legend()

    axs[2,0].plot(cell.Qconc_data, label="Qconc")
    axs[2,0].set_title('Flujo de concentrado')
    axs[2,0].legend()

    axs[2,1].plot(cell.valve_aperture, label="u1")
    axs[2,1].set_title("Apertura Válvula")
    axs[2,1].legend()

    axs[2,2].plot(cell.air_flux, label="u2")
    axs[2,2].set_title("Flujo de aire")
    axs[2,2].legend()
    plt.tight_layout()
    plt.legend()
    plt.show()
    

def main(): 
    # Number of inputs: 
    Ninputs = 4
    mv_time_span = []
    for i in range(0, Ninputs):
        u1 = float(input("Valor de apertura de válvula: "))
        u2 = float(input("Valor de flujo de aire: "))
        mv_time_span.append([u1, u2])

    # Cell instantiation: 
    cell_id = 1
    i = 2
    K = 5
    tspan = (0, 4e4)
    u_time_interval = int(tspan[-1] / 4)
    cell = Cell(cell_id=cell_id, i=i, num_classes=K, 
                time_span=tspan, u_time_span=mv_time_span, 
                u_time_interval=u_time_interval)
    
    start_time = time.time()
    cell.update_model()
    c = cell.iter_counter
    print(f"Iterations: {c}")
    end_time = time.time()
    print(f"Code execution time: {round(end_time - start_time)} [s]")

    plot_control_variables(cell=cell)


if __name__ == "__main__": 
    main()




