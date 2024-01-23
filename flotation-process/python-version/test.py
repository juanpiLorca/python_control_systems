import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# Parametros y constantes globales
i = 2
a = -6.06e-4
b = 0.9182
c = 0
n = 1
C = 6.38e-4
aSAT = 0.01
bSAT = 0.4927
cSAT = 10.48
dSAT = 0.2588
aSF = 50
bSF = 0.5
klambda = 6.815
g = 9.81
CPB = 50
Pe = 0.15
hT = 0.48
kv = 1e-3
f = 0.8
K = 5
Llip = 1.51
hover = 0.15
freq_bubble_size = np.array([500, 1300, 1700, 2000, 1200])
db_k_pulp = np.array([2e-4, 4e-4, 6e-4, 8e-4, 1e-3])
tspan = [0, 30]
epsis = np.array([0.18, 0.16, 0.14, 0.12, 0.1])
x0 = np.concatenate((np.zeros(i), [0.13], epsis / (1 - np.sum(epsis)), [1e-6]))
y0 = np.concatenate(([0], np.zeros(i), np.zeros(i), [0], [0]))

Cfeed = np.array([0.2,0.8])
Qfeed = 1.67e-4
P = np.array([2.8e-4,3e-4])
dp = np.array([75e-6,73e-6])
Vcell = 0.087
Acell = 0.18


phi = 0.2
density_sol = 2500
density_water = 1000
viscosity_water = 1e-3
density_pulp = phi * density_sol + (1 - phi) * density_water
viscosity_pulp = viscosity_water * np.exp(2.5 * phi / (1 - 0.609 * phi))



# Valores iniciales
yant = np.zeros(7)
u = [0.3, 1e-3]  # [0.25, 7e-4]


# Funcion que calcula las derivadas del modelo
def odefun(t,x):
    #global i, u, a, b, c, n, C, aSAT, bSAT, cSAT, dSAT, aSF, bSF, kv, Cfeed, Qfeed, Vcell, Acell, P, density_sol, density_water, phi, density_pulp, viscosity_pulp, klambda, g, CPB, Pe, hT, f, K, freq_bubble_size, db_k_pulp, dp, Llip, hover
    #global yant

    m = x[:i]
    h0 = x[i]
    gh = x[i+1:i+1+K]
    dbf = x[i+1+K]

    s = model_equations(x)
    gh_tot = s[0]
    hp = s[1]
    Qair_k = s[2:7]
    vg_k_pulp_out = s[7:12]
    mfeed = s[12:14]
    mtail = s[14:16]
    mTF = s[16:18]
    mENT = s[18:20]
    Qconc = s[20]
    Qtail = s[21]
    alpha_sat = s[22]
    not_solve = s[23]

    # Derivadas del estado

    DdbfDt = C * dbf**(1 - n)

    if not not_solve:
        Qpulp_out = Qtail + Qconc

        Dh0Dt = (1 / Acell) * (Qfeed - Qpulp_out)
        DghDt = ((1 + gh_tot) / (Acell * hp)) * (Qair_k - Acell * vg_k_pulp_out * (gh / (1 + gh_tot)) - (Qfeed - Qpulp_out) * gh)
        DghDt = DghDt.reshape((K,))
        DmDt = mfeed - mtail - mTF - mENT

        #y = np.concatenate(([alpha_sat],mTF,mENT,[Qtail],[Qconc]))

    else:
        Dh0Dt = 0
        DghDt = np.zeros(K)
        DmDt = np.zeros(i)
        #y = yant

    #yant = y
    #Y.append(y)

    return np.concatenate((DmDt, [Dh0Dt], DghDt, [DdbfDt]))


# Funcion que define las ecuaciones del modelo
def model_equations(x):

    m = x[:i]
    h0 = x[i]
    gh = x[i+1:i+1+K]
    dbf = x[i+1+K]

    gh_tot = np.sum(gh)
    epsi0 = gh / (1 + gh_tot)

    uv = u[0]  # valvula flujo de colas
    Qair_in = u[1]  # flujo de alimentacion de aire

    gh_tot = np.sum(gh)
    epsi0 = gh / (1 + gh_tot)
    epsi0_tot = np.sum(epsi0)
    hp = h0 / (1 - epsi0_tot)

    hf = hT - hp
    vg_k_pulp_out = (g * density_pulp * db_k_pulp**2) / (18 * viscosity_pulp * (1 - epsi0)**1.39)
    db_int = np.sum(vg_k_pulp_out * epsi0) / np.sum(vg_k_pulp_out * epsi0 / db_k_pulp)
    jg = Qair_in / Acell
    vb = a + b * jg + c * jg**2
    K1 = density_pulp * g / (3 * CPB * viscosity_pulp)
    Qtail = kv * uv * np.sqrt(h0)
    vg_total_out_pulp = np.sum(vg_k_pulp_out * epsi0)

    psi_k = freq_bubble_size / np.sum(freq_bubble_size)
    Qair_k = Qair_in * psi_k

    # Resolvemos el loop algebraico
    def equations(vars):
        vg_, alpha = vars
        alpha_ = (vg_ - vb) / vg_
        alpha_sat = aSAT + (bSAT / (1 + np.exp(-cSAT * (alpha_ - dSAT))))
        tau_f = hf / vg_
        dbf_out = (n * C * tau_f + db_int**n)**(1/n)
        lambda_out = klambda / dbf_out**2
        Sig = 1 / (1 + np.exp(-aSF * (alpha - bSF)))
        Qconc = ((1 - Sig) * (Acell * vg_**2 * lambda_out) * (1 - alpha_sat) * alpha_sat / K1) + (Sig * Acell * vg_**2 * lambda_out / 4 * K1)
        DghDt = ((gh_tot + 1) / (Acell * hp)) * (Qair_k - Acell * vg_k_pulp_out * (gh / (1 + gh_tot)) - (Qfeed - Qtail - Qconc) * gh)
        epsi = ((1 - Sig) * vg_ * (1 - alpha_sat) * lambda_out / K1) + (Sig * vg_ * lambda_out / (2 * K1))
        eqn1 = vg_ - (h0 * np.sum(DghDt) - ((gh_tot + 1) * (Qfeed - Qtail - Qconc)) / Acell) - vg_total_out_pulp
        eqn2 = alpha - (Qconc / (epsi * Qair_in))
        return [eqn1, eqn2]

    # Initial guess for the solution
    initial_guess = [1e-5, 1e-5]
    solution = fsolve(equations, initial_guess)

    if not np.any(np.isnan(solution)):
        not_solve = 0
        vg_, alpha = solution
        alpha_ = (vg_ - vb) / vg_
        alpha_sat = aSAT + (bSAT / (1 + np.exp(-cSAT * (alpha_ - dSAT))))
        tau_f = hf / vg_
        dbf_out = (n * C * tau_f + db_int**n)**(1/n)
        lambda_out = klambda / dbf_out**2
        Sig = 1 / (1 + np.exp(-aSF * (alpha - bSF)))
        Qconc = ((1 - Sig) * (Acell * vg_**2 * lambda_out) * (1 - alpha_sat) * alpha_sat / K1) + (Sig * Acell * vg_**2 * lambda_out / 4 * K1)
    else:
        not_solve = 1

    if not not_solve:
        Daxial = (jg**1.5) / (np.sqrt(K1 * (np.sqrt(3) - np.pi / 2)) * Pe)
        Sb = 6 * vg_ / db_int
        vset = g * (density_sol - density_water) * dp**2 * (1 - phi)**4.65 / (18 * viscosity_pulp * 3)

        RF = ((1 - Sig) * (((alpha_sat * (1 - alpha_sat) * vg_) / vset)**(f / 2)) * ((db_int / dbf)**f)) + Sig * ((vg_ / vset)**(f / 2)) * ((db_int / dbf)**f)
        ENT = ((1 - Sig) * np.exp(-(vset**1.5) * hf / (Daxial * np.sqrt(vg_ * (1 - alpha_sat))))) + Sig * np.exp(-2 * (vset**1.5) * hf / (Daxial * np.sqrt(vg_)))

        k = P * Sb

    
    Vk_gas = (epsi0 / (1 - epsi0)) * h0 * Acell
    Vgas = np.sum(Vk_gas)
    Vpulp = h0 * Acell + Vgas
    Ctail = m / Vpulp


    # Masas 
    mfeed = Cfeed * Qfeed
    mtail = Ctail * Qtail
    mTF = Vcell * k * RF * Ctail
    mENT = Qconc * ENT * Ctail

    return np.concatenate(([gh_tot], [hp], Qair_k, vg_k_pulp_out, mfeed, mtail, mTF, mENT, [Qconc] , [Qtail], [alpha_sat], [not_solve]))


# Resolvemos las ecuaciones diferenciales
#solution = solve_ivp(odefun, tspan, x0, method='RK45',events=event)
solution = solve_ivp(odefun, tspan, x0, method='RK45')

model_vars = []
for k in range(len(solution.y[0,:])):    
    model_vars.append(model_equations(solution.y[:,k]))
# model_vars[0][1] es la solucion de indice 0 y el elemento 1 de dicha solucion

hp = []
for k in range(len(model_vars)):
    hp.append(model_vars[k][1])

mTF = []
for k in range(len(model_vars)):
    mTF.append(model_vars[k][16])

mENT = []
for k in range(len(model_vars)):
    mENT.append(model_vars[k][18])

alpha_sat = []
for k in range(len(model_vars)):
    alpha_sat.append(model_vars[k][22])

Qtail = []
for k in range(len(model_vars)):
    Qtail.append(model_vars[k][21])

Qconc = []
for k in range(len(model_vars)):
    Qconc.append(model_vars[k][20])



#print(solution.t_events)
#print(solution.t_events[0])
#print(solution.y_events)
#print(Y)


# ----------------------------------------------------- Graficos

fig, ax = plt.subplots(3,2)

ax[0,0].plot(solution.t, solution.y[0,:]) # estado 0
ax[0,0].plot(solution.t, solution.y[1,:]) # estado 1
ax[0,0].set_title('Masa total de cada mineral')

ax[0,1].plot(solution.t, mTF)
ax[0,1].plot(solution.t, mENT)
ax[0,1].set_title('Mtf y Ment interes economico')

ax[1,0].plot(solution.t, hp)
ax[1,0].set_title('Altura de la pulpa en la celda')

ax[1,1].plot(solution.t, alpha_sat)
ax[1,1].set_title('Recuperacion de aire')

ax[2,0].plot(solution.t, Qtail)
ax[2,0].set_title('Flujo de cola')

ax[2,1].plot(solution.t, Qconc)
ax[2,1].set_title('Flujo de concentrado')

plt.tight_layout() 
plt.show()    