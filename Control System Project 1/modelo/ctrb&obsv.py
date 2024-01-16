from control.matlab import *
import numpy as np 

# ---- Controlabilidad ----

# Definimos la matriz de transición de estados A evaluada en el punto de eq. genérico: 
A = np.matrix("0 0 1 0; 0 0 0 1; 0 0 0 0;  0 0 0 -0.20833")

# Definimos la matriz de controles B evaluada en el punto de eq. genérico: 
B = np.matrix("0 0; 0 0; 0.000022727 0; 0 0.000020833")

# Calculamos la matriz de controlabilidad: 
C = ctrb(A, B)

# Verificamos que el rango de tal matriz sea igual a la cantidad de variables de estado: rank(C) = 4
rank_c = np.linalg.matrix_rank(C)
print(f"El rango de la matriz de controlabilidad es: {rank_c}")
# Por lo que el rango de tal matriz es igual a la cantidad de variables de estado, entonces el sistema es controlable. 


# ---- Observabilidad ----

# Mediciones x_dot y theta_dot por parte del sensor: 
C1 = np.matrix(" 0 0 1 0; 0 0 0 1")

# Calculamos la matriz de observabilidad: 
O1 = obsv(A, C1)
print(O1)

# Verificamos que el rango de tal matriz sea igual a la cantidad de variables de estado: rank(C) = 4
rank_o1 = np.linalg.matrix_rank(O1)
print(f"El rango de la matriz de observabilidad para las mediciones x_dot y theta_dot es: {rank_o1}")
# Por lo tanto no es observable para las mediciones de velocidad del carro y velocidad angular cañón: rank(O1) = 2 != 4


# Mediciones x y theta por parte del sensor: 
C2 = np.matrix(" 1 0 0 0; 0 1 0 0")

# Calculamos la matriz de observabilidad: 
O2 = obsv(A, C2)
print(O2)

# Verificamos que el rango de tal matriz sea igual a la cantidad de variables de estado: rank(C) = 4
rank_o2 = np.linalg.matrix_rank(O2)
print(f"El rango de la matriz de observabilidad para las mediciones x_dot y theta_dot es: {rank_o2}")
# Por lo tanto si es observable para las mediciones de posición del carro y posición angular cañón: rank(O2) = 4 = 4


