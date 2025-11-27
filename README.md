# Proyecto: Simulación del Método de Newton–Raphson

Integrantes:


--Samuel Arenas Valencia 2341928 

--Nicolás David Córdoba 2343576

--Daniel Andrade Reyes 2343792

Este proyecto implementa una simulación numérica del método de Newton–Raphson para resolver sistemas de ecuaciones no lineales. A continuación encontrarás una descripción clara y sencilla de los archivos principales y de cómo se conecta todo dentro del proyecto.

## Archivos del proyecto

### 1. `vector.py`
Aquí se encuentra la clase `NewtonRaphson`, responsable de toda la parte matemática del método. Entre sus funciones principales están:
- Definir el vector inicial (`vecInicial`) con las condiciones de frontera y los valores de los nodos.
- Calcular el vector de residuales (`cal_function`), es decir, evaluar las ecuaciones del sistema en cada iteración.
- Construir la matriz Jacobiana (`cal_jacobiano`) para actualizar la solución.
- Actualizar el vector de solución usando gradiente conjugado o la inversa de la Jacobiana (`newVector`, `newVectorInversa`).
- Analizar propiedades de la Jacobiana: simetría, diagonal dominante, rango completo, etc.
- Interpolar con splines cúbicos mediante `construir_spline_natural` para obtener curvas más suaves del campo de velocidades.

### 2. `plotVec.py`
Este archivo contiene la clase `Plot`, encargada de generar las visualizaciones de la simulación. Permite:
- Mostrar mapas de calor de la solución (`showPlot`).
- Crear animaciones de la evolución del sistema (`animateHistoryLoop`).
- Generar gráficos suavizados mediante interpolación bicúbica (`showPlotSpline`).
- Ver detalles específicos del dominio con funciones como `showPlotDetail` o `showPlotIndices`.

### 3. `interfaz.py`
Aquí se implementa `SimulationViewer`, la interfaz gráfica del proyecto. Su función es permitirte explorar los resultados de la simulación de forma visual e interactiva:
- Modo de “diapositivas” para recorrer distintos tipos de gráficos.
- Integración con Matplotlib para mostrar mapas de calor, animaciones y gráficos interpolados dentro de la ventana.
- Botones de navegación para analizar las vistas a tu ritmo.
- Adaptadores visuales como `view_heatmap`, `view_spline_static`, `view_anim_history` y `view_anim_spline`.

### 4. `main.py`
Es el punto de entrada del proyecto y coordina el flujo de la simulación:
- Configura la malla (`get_scaled_config`) ajustando escalas y condiciones iniciales.
- Ejecuta las iteraciones de Newton-Raphson (`iteration`) usando gradiente conjugado o la inversa de la Jacobiana.
- Genera distintas configuraciones de simulación para comparar resultados.
- Lanza la interfaz gráfica con `SimulationViewer` para visualizar todo de forma cómoda.

## Cómo ejecutar el proyecto

1. Instala las dependencias necesarias (`numpy`, `matplotlib`, `scipy`).
2. Ejecuta el archivo principal:

python main.py