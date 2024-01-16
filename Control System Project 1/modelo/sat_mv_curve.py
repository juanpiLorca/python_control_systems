import numpy as np 
import matplotlib.pyplot as plt 

# Chart's engine saturation curve [N]:
x1 = np.linspace(-50000, 50000, 100000)
y1 = np.zeros_like(x1)

y1[:20000] = -30000
y1[20000:80000] = x1[20000:80000]
y1[80000:] = 30000

fig1, ax1 = plt.subplots()
ax1.plot(x1, y1)
plt.title("Curva de saturación fuerza de motor carro [N]")
plt.grid()
plt.show()

# Cannon's engine saturation curve [Nm]:
x2 = np.linspace(-300000, 300000, 600000)
y2 = np.zeros_like(x2)

y2[:64560] = -235440
y2[64560:535440] = x2[64560:535440]
y2[535440:] = 235440

fig2, ax2 = plt.subplots()
ax2.plot(x2, y2)
plt.title("Curva de saturación torque de elevación [Nm]")
plt.grid()
plt.show()