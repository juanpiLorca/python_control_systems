import matplotlib.pyplot as plt 
import matplotlib
import pandas as pd

matplotlib.use('TkAgg')  # Puedes probar 'Qt5Agg' u otros backends interactivos si 'TkAgg' no funciona

from matplotlib.animation import FuncAnimation, Animation

t_vals = []
y_vals = []

def plot_real_time(i): 
    data = pd.read_csv("data.csv")
    t = data["time"]
    y = data["y_value"]
    
    plt.cla()
    plt.plot(t, y, label="Model Data")
    plt.legend(loc="upper left")
    plt.tight_layout()

realTime = FuncAnimation(plt.gcf(), plot_real_time, interval=500, cache_frame_data=False)


plt.tight_layout()
plt.show()
