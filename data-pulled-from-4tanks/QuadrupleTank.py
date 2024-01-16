import numpy as np 
import pandas as pd 
import torch
import os 

import matplotlib.pyplot as plt


class QuadrupleTank(): 
    """
    Linearized around x0 state-space quadruple-tanks system
    args: 
        >>> x0: operating point 
        >>> Ts: sampling time
        >>> Kr: reference gain
    """
    
    def __init__(self, x0: list, Ts: float, Kr: float):
        self.x0 = np.array(x0)
        self.x = np.array(x0)
        self.u = np.zeros((2,))
        self.ref = np.zeros((2,))
        self.Ts = Ts
        self.current_time = 0

        # State-space system
        self.A = np.array([ [-0.0089, 0,          0.0089,     0       ],
                            [0,       -0.0062,    0,          0.0101  ],
                            [0,       0,          -0.0089,    0       ],
                            [0,       0,          0,          -0.0062 ] ])  
        self.B = np.array([ [0.0833,  0       ],
                            [0,       0.0624  ],
                            [0,       0.0476  ],
                            [0.0312,  0       ] ])   
        self.C = np.array([ [1, 0,  0,  0],
                            [0, 1,  0,  0]    ])  
        self.D = np.zeros((2, 2)) 

        # Use LQR to design the state feedback gain matrix K
        self.K = np.array([ [0.7844,    0.1129,   -0.0768,    0.5117],
                            [0.0557,    0.7388,    0.5409,   -0.0397] ])
        self.Kr = Kr
        # Closed-loop dynamics (A - BK)
        self.CL_dyn = self.A - self.B @ self.K

        # Storing variables
        self.x_states = []
        self.u_data = []

    # Open loop dynamics: input is u (2x1) and output is x (4x1)
    def open_loop(self, u):
        # Calculate the deviation from equilibrium
        dx = self.A @ self.x + self.B @ u

        # Update of the state
        self.x = self.x + self.Ts * dx

        # Increment time
        self.current_time += self.Ts

    # Closed loop dynamics: states and inputs are saved iteratively in x and u, respectively
    def closed_loop(self):
        # Calculate control input u = -Kx
        self.u = self.Kr * self.ref - self.K @ self.x

        # Update state: x_dot = (A - B.K).x + B.r * Kr
        self.x = self.x + self.Ts * (self.CL_dyn @ self.x + self.Kr * self.B @ self.ref)

        # Increment time
        self.current_time += self.Ts

    def store(self): 
        # Storing data: 
        self.x_states.append(self.x)
        self.u_data.append(self.u)



def gen_inputs(n_seqs, seqlen): 
    "Inputs: voltage values"
    pass

def gen_references(n_seqs, seqlen): 
    "References: water height reference in each tank"
    ref1_seqs = []
    ref2_seqs = []

    for i in range(n_seqs): 
        ref1 = max(np.clip(np.random.normal(size=(1,))[0], 0, 1), 0.05)
        ref2 = max(np.clip(np.random.normal(size=(1,))[0], 0, 1), 0.05)
        ref1_seqs.append([ref1 for j in range(seqlen)])
        ref2_seqs.append([ref2 for j in range(seqlen)])
    
    ref1_seqs = np.concatenate(ref1_seqs).reshape(-1, 1)
    ref2_seqs = np.concatenate(ref2_seqs).reshape(-1, 1)
    return ref1_seqs, ref2_seqs




if __name__ == "__main__": 
    ref1, ref2 = gen_references(n_seqs=3000, seqlen=60)
    N = len(ref1)
    ref1, ref2 = np.array(ref1), np.array(ref2)
    refs = np.concatenate([ref1, ref2], axis=1)
    ref1 *= 100
    ref2 *= 100
    x0 = [40, 40, 40, 40]
    Ts = 0.01
    Kr = 1.0
    system = QuadrupleTank(x0=x0, Ts=Ts, Kr=Kr)

    for i in range(0, N): 
        r1 = ref1[i]
        r2 = ref2[i]
        ref = np.array([r1, r2])
        system.ref = ref
        system.closed_loop()
        system.store()

    time_series = np.array(system.x_states)
    time_series = time_series[:, :, 0]
    t = np.linspace(0, system.current_time, time_series.shape[0])
    print(f"Time series shape (plant states): {time_series.shape}")

    fig, ax = plt.subplots(nrows=4, ncols=1)
    ax[0].plot(t, time_series[:, 0], label="x_0")
    ax[0].grid()
    ax[0].legend()

    ax[1].plot(t, ref1, "--", label="ref1", color="C2")
    ax[1].grid()
    ax[1].legend()

    ax[2].plot(t, time_series[:, 1], label="x_1")
    ax[2].grid()
    ax[2].legend()

    ax[3].plot(t, ref2, "--", label="ref2", color="C2")
    ax[3].grid()
    ax[3].legend()

    # ax[4].plot(t, time_series[:, 2], label="x_2")
    # ax[4].grid()
    # ax[4].legend

    # ax[5].plot(t, time_series[:, 3], label="x_3")
    # ax[5].grid()
    # ax[5].legend()

    plt.tight_layout()
    plt.legend()
    plt.show()

    # Data storing:
    data = np.concatenate([refs, time_series], axis=1)
    print(f"Data to store shape: {data.shape}")

    # mkdir:
    if not os.path.exists('No_Noised_Inputs/'):
        os.makedirs('No_Noised_Inputs/')

    # NaN values check:
    has_nan = np.isnan(data).any()
    if has_nan:
        print("The matrix contains NaN values.")
        data = np.array([[0 if np.isnan(element) else element for element in row] for row in data])

    has_nan = np.isnan(np.array(data)).any()
    if has_nan:
        print("The matrix STILL contains NaN values.")
    else:
        print("NaN values were deleted.")

    np.save('No_Noised_Inputs/data_fcn_clean.npy', data)
    #np.save('data_noise.npy', np.concatenate([inputs_noise, series_noise], axis=1))
    torch.save(data, 'No_Noised_Inputs/data_clean.pkl')
    print(f'Saved \'data\'')



