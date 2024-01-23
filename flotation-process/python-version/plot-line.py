import matplotlib.pyplot as plt 
import pandas as pd 

data = pd.read_csv("data_line.csv")
fieldnames = data.columns.tolist()

t = data["time"]
N = len(fieldnames) - 1

fig, axs = plt.subplots(nrows=N, ncols=1)
for i in range(N): 
    x = data[fieldnames[i+1]]
    axs[i].plot(t, x, label=fieldnames[i+1])
    axs[i].set_title(f"Variable: {fieldnames[i+1]}")
    axs[i].legend()
    axs[i].grid(True)

plt.tight_layout()
plt.legend()
plt.show()

