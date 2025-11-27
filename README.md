# Proyecto: Simulación Numérica de Flujo Incompresible utilizando las Ecuaciones de Navier–Stokes

Integrantes:

--Samuel Arenas Valencia 2341928 

--Nicolás David Córdoba 2343576

--Daniel Andrade Reyes 2343792


Este proyecto implementa una simulación numérica basada en el método de Newton-Raphson para resolver sistemas de ecuaciones no lineales. A continuación, se describe la función de cada archivo en el proyecto:

---

## Archivos del Proyecto

### 1. `vector.py`
Este archivo contiene la clase `NewtonRaphson`, que implementa toda la lógica matemática y numérica del método de Newton-Raphson. Sus principales responsabilidades incluyen:
- **Inicialización del vector de simulación (`vecInicial`)**: Establece las condiciones iniciales del sistema, como los valores en los nodos y las condiciones de frontera.
- **Cálculo del vector de residuales (`cal_function`)**: Evalúa las ecuaciones discretizadas en cada nodo para determinar el error actual.
- **Cálculo de la matriz Jacobiana (`cal_jacobiano`)**: Construye la matriz Jacobiana del sistema para cada iteración.
- **Métodos de actualización del vector (`newVector`, `newVectorInversa`)**: Calcula el nuevo vector de solución utilizando el gradiente conjugado o la inversa de la Jacobiana.
- **Propiedades de la matriz Jacobiana**: Incluye métodos para verificar si la matriz es simétrica, diagonalmente dominante, de rango completo, entre otros.
- **Interpolación con splines cúbicos (`construir_spline_natural`)**: Implementa un spline cúbico natural para interpolar los valores del vector de velocidades.

---

### 2. `plotVec.py`
Este archivo contiene la clase `Plot`, que se encarga de la visualización de los resultados de la simulación. Sus principales responsabilidades incluyen:
- **Mapas de calor (`showPlot`)**: Genera gráficos de los valores del vector de simulación en forma de mapa de calor.
- **Animaciones de la simulación (`animateHistoryLoop`)**: Muestra la evolución del sistema a lo largo de las iteraciones.
- **Interpolación y suavizado (`showPlotSpline`)**: Aplica interpolación bicúbica para generar gráficos suavizados de los resultados.
- **Visualización detallada (`showPlotDetail`, `showPlotIndices`)**: Muestra los valores exactos o los índices de los nodos en la malla para facilitar el análisis.

---

### 3. `interfaz.py`
Este archivo implementa la clase `SimulationViewer`, que proporciona una interfaz gráfica interactiva para visualizar los resultados de la simulación. Sus principales características incluyen:
- **Visor de diapositivas**: Permite navegar entre diferentes vistas de la simulación, como mapas de calor, animaciones y gráficos interpolados.
- **Integración con Matplotlib**: Utiliza Matplotlib para incrustar gráficos y animaciones en la interfaz.
- **Controles de navegación**: Incluye botones para avanzar o retroceder entre las vistas.
- **Adaptadores de visualización**: Define funciones como `view_heatmap`, `view_spline_static`, `view_anim_history` y `view_anim_spline` para conectar los datos de la simulación con los gráficos.

---

### 4. `main.py`
Este archivo es el punto de entrada del proyecto. Contiene la lógica principal para configurar y ejecutar las simulaciones. Sus principales responsabilidades incluyen:
- **Configuración de la malla (`get_scaled_config`)**: Escala el tamaño de la malla y ajusta las condiciones iniciales según un factor de escala.
- **Ejecución del método de Newton-Raphson (`iteration`)**: Realiza las iteraciones del método utilizando el gradiente conjugado o la inversa de la Jacobiana.
- **Instanciación de simulaciones**: Crea múltiples configuraciones de simulación con diferentes factores de escala y condiciones de frontera.
- **Visualización de resultados**: Utiliza la clase `SimulationViewer` para mostrar los resultados de las simulaciones en una interfaz gráfica.

Detalles y uso de `get_scaled_config(factor)`:
- Firma: `get_scaled_config(factor: float) -> (Nx_new, Ny_new, h_new, blocks_new)`
- Qué hace:
  - Parte de una configuración base: `Nx_base = 81`, `Ny_base = 9`, `h_base = 5.0`, `blocks_base = [[10, 20, 0, 2], [58, 80, 5, 8]]`.
  - Calcula Nx_new = int(Nx_base * factor) y Ny_new = int(Ny_base * factor).
  - Ajusta el paso espacial a `h_new = h_base / factor` (factor > 1 → mayor resolución).
  - Escala las coordenadas de los bloques internamente y corrige las coordenadas superiores para que no excedan los límites de la nueva malla.
  - Devuelve `(Nx_new, Ny_new, h_new, blocks_new)`.

- Ejemplos de uso:
  - Factor 1 (malla base): `nx, ny, h, blocks = get_scaled_config(1)`
  - Factor 3 (malla 3× en cada dirección): `nx, ny, h, blocks = get_scaled_config(3)`
  - Factor 0.5 (malla reducida): `nx, ny, h, blocks = get_scaled_config(0.5)`

- Consideraciones prácticas:
  - Aumentar `factor` aumenta Nx y Ny y reduce `h` → mejora precisión pero incrementa memoria y tiempo de cómputo.
  - Para pruebas rápidas use factor pequeño o reduzca `Ny_base`/`Nx_base`.
  - El escalado de `blocks` ajusta las coordenadas de los bloques internamente; si modificas `blocks_base`, asegúrate de mantenerlos consistentes con la malla.
  - Puedes cambiar `h_base` directamente dentro de `get_scaled_config` si quieres un valor por defecto distinto al del archivo.

Parámetros configurables en `main.py` y su influencia:
- En instanciación de `NewtonRaphson(...)`:
  - `V0`: valor máximo aparente para visualización.
  - `v_ij`: valor de velocidad inicial dentro del dominio.
  - `v0ParedSuperior`: condición de frontera en la pared superior.
  - `h`: paso espacial (retornado por `get_scaled_config`).
  - `blocks`: bloques escalados para condiciones internas.
- En `iteration()`:
  - `max_iterations`: máximo de iteraciones.
  - `epsilon`: tolerancia para el residuo (criterio de parada).
  - `delta`: tolerancia para el cambio entre iteraciones.
  - `saveHistory`: si `True`, guarda el historial para animaciones/visualización.
  - `callNewVector`: indica si se aplica gradiente conjugado o inversa de jacobiano.

Sugerencias para experimentación:
- Comparar convergencia y rendimiento:
  - Cambia `FACTORA` y `FACTORB` para ejecutar pruebas A/B con distintas resoluciones.
  - Mantén una simulación de referencia con `factor=1` y otra con `factor=3` para observar diferencias.
- Prueba condiciones de frontera y velocidad inicial:
  - Modifica `v0ParedSuperior`, `V0` y `v_ij` en las instancias `NewtonRaphson(...)` para distintos escenarios.
- Visualización y animaciones:
  - Activa `saveHistory=True` en `iteration()` para generar animaciones de la evolución.
  - Usa `SimulationViewer` para alternar entre vistas (mapa de calor, spline, animaciones).