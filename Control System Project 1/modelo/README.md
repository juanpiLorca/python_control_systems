# modelo

modelo es la carpeta utilizada para la simulación del proyecto de control automático: Cañón Contra Incendios

## Files:

La carpeta modelo presenta dos tipos de archivos: 
1.- Utilizados en la simulación: 
    * model.py
    * fire.py 
    * projectile.py
    * sim.py
    * plot_data.py
2.- Utilizados en los cálculos matemáticos realizados para los distintos puntos de enunciado.
    * ctrb&obsv.py
    * saturation_curve.py

## Modo de uso: 

Para la ejecución del modelo se debe correr el script sim.py que contiene la clase ModelDisplayPixel() que realiza la simulación del proyecto. Todos los cálculos hechos para la posición horizontal y el ángulo de elevación son hechos en el script model.py que contiene la clase Model(). Para los cálculos de la trayectoria de la bombarda, se utiliza el script projectile.py que contiene la clase Projectile(). Por último, para la generación de fuegos aleatorios, está el script fire.py que contien la clase Fire(). Todos los scripts encargados de hacer cálculos matemáticos para la simulación envían información a la simulación principal en sim.py.

## ¡Indicaciones importantes!

==> Storing datapoints: 
    El modelo logra guardar puntos relevantes de la simulación a través del tiempo y graficarlos en el script plot_data.py. El modo de guardado es AUTOMÁTICO a los 100 segundos de simulación. Si se desea guardar los datos en menos tiempo se debe cambiar la variable k_data_samples ubicada en la línea 53 de código al valor que usted determine conveniente (ej: k_data_samples = 2000 que son 40 segundos de simulación hasta el primer guardado) 

==> Punto graficado al hacer click: 
    El sistema recibe como referencia el punto graficado en la pantalla (PUNTO BLANCO), por lo que regula el movimiento del carro y cañón de acuerdo a la posición fijada por el usuario. NOTAR: que el ángulo referencial es referido al centro de masa inicial del carro menos 2000 kilómetros o 100 pixeles, que es dónde se pararía el carro a disparar.

==> Curvas de saturación variables manipuladas: 
    Para la sección del controlador: x_controller() y theta_controller() pueden comentar la "saturation curve" para ver como responde el sistema sin márgenes que limiten estas variables.

==> Debouncer: 
    Mucho ojo con cómo se apreta la tecla A. Es necesario verificar el mensaje en la consola que dice que el modo está en automático. 

==> Perturbaciones y ruido en los sensores para la medición de theta: 
    Se puede descomentar el "# + np.random.normal(0,1)" para ver los efectos de estas perturbaciones en la respuesta del sistema para theta.


