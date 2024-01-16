import numpy as np
import matplotlib.pyplot as plt 

Ts = 0.1
def IQE(error_array, Ts, ref, variable): 
    print(f"Índice de desempeño IQE para {variable} con referencia {ref}")
    print(f"IQE = {np.sum(np.square(error_array*Ts))}") 

print("Loading data...")
data = np.load("modelo/simdata.npy")
print(f"data shape: {(data.shape[0],data.shape[1])}")

t = data[0]
x = data[1]
theta = data[2]
x_ref = data[3]
theta_ref = data[4]
mv_x = data[5]
mv_c = data[6]

# IQE x-variable:
IQE((x-x_ref), Ts, x_ref[len(x_ref)-1], "x")

# IQE theta-variable:
IQE((theta-theta_ref), Ts, theta_ref[len(theta_ref)-1], "theta")

print("Drawing plot 1...")
fig1 = plt.figure()
plt.title('State variables: horizontal cart position, Reference')
plt.plot(t, x, t, x_ref)
plt.grid()
plt.xlabel('t [s]')
plt.ylabel('Position [m]')

print("Drawing plot 2...")
fig2 = plt.figure()
plt.title('State variables: cannon elevation, Reference')
plt.plot(t, theta, t, theta_ref)
plt.grid()
plt.xlabel('t [s]')
plt.ylabel('Angle [rad]')

print("Drawing plot 3...")
fig3 = plt.figure()
plt.title('State error: x-postion')
plt.plot(t, x_ref-x)
plt.grid()
plt.xlabel('t [s]')
plt.ylabel('error [m], ref [m]')

print("Drawing plot 4...")
fig4 = plt.figure()
plt.title('State error: theta-elevation')
plt.plot(t, theta_ref-theta)
plt.grid()
plt.xlabel('t [s]')
plt.ylabel('error [rad], ref [rad]')

print("Drawing plot 5...")
fig5 = plt.figure()
plt.title('Manipulate Variable: engine force')
plt.plot(t, mv_x)
plt.grid()
plt.xlabel('t [s]')
plt.ylabel('engine [N]')

print("Drawing plot 6...")
fig6 = plt.figure()
plt.title('Manipulate Variable: cannon engine torque')
plt.plot(t, mv_c)
plt.grid()
plt.xlabel('t [s]')
plt.ylabel('engine torque [Nm]')

plt.show()

