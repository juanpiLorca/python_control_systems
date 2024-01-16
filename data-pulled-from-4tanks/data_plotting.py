import numpy as np
import matplotlib.pyplot as plt 
import pickle
from data_processor import generate_noise


#pkl_file_path = 'No_Noised_Inputs/data_dict_processor_{}.pkl'
data_clean_path = 'No_Noised_Inputs/data_fcn_clean.npy'
data_clean = np.load(data_clean_path)

noisy_data = generate_noise(data_clean=data_clean, plot=True)

hTanks = data_clean[:, 2:4]
noisy_hTanks = noisy_data[:, 2:4]

fig, ax = plt.subplots(nrows=2, ncols=1)
ax[0].plot(noisy_hTanks[:, 0], label="h1 noisy")
ax[0].plot(hTanks[:, 0], label="h1 clean")
ax[0].grid()
ax[0].legend()
ax[1].plot(noisy_hTanks[:, 1], label="h2 noisy")
ax[1].plot(hTanks[:, 1], label="h2 clean")
ax[1].grid()
ax[1].legend()
plt.tight_layout()
plt.legend()
plt.show()