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

---

## Cómo Ejecutar el Proyecto

### Requisitos Previos

Asegúrate de tener instalado Python 3.x y las siguientes bibliotecas:

```bash
pip install numpy matplotlib scipy
```

### Ejecutar la Simulación

Para correr el proyecto, simplemente ejecuta el archivo principal:

```bash
python main.py
```

Esto iniciará las simulaciones configuradas y abrirá una interfaz gráfica interactiva (`SimulationViewer`) donde podrás navegar entre diferentes visualizaciones de los resultados.

---

## Configuración de Parámetros

### 1. Cambiar el Factor de Escalado del Sistema

En `main.py`, localiza las siguientes líneas (aproximadamente línea 125):

```python
FACTORA = 1
FACTORB = 3   # Hacer el factor más pequeño hace que la prueba B corra más rápido
```

- **`FACTORA`**: Factor de escala para las simulaciones del grupo A (por defecto: 1)
- **`FACTORB`**: Factor de escala para las simulaciones del grupo B (por defecto: 3)

**Impacto del factor:**

- `factor = 1`: Malla base (81×9 nodos, rápida)
- `factor = 3`: Malla ampliada (243×27 nodos, más precisa pero más lenta)
- `factor = 0.5`: Malla reducida (40×4 nodos, muy rápida para pruebas)

**Recomendación:** Para pruebas rápidas, usa factores pequeños (0.5, 1). Para resultados detallados, usa factores mayores (2, 3, 5).

### 2. Modificar las Condiciones Iniciales del Sistema

Dentro de la función `get_scaled_config(factor)` en `main.py` (líneas 8-45), puedes modificar:

```python
# 1. Configuración Base (Original)
Nx_base = 81        # Número de nodos en dirección X
Ny_base = 9         # Número de nodos en dirección Y
h_base = 5.0        # Paso espacial base
blocks_base = [[10, 20, 0, 2], [58, 80, 5, 8]]  # Bloques internos [x1, x2, y1, y2]
```

**Parámetros modificables:**

- **`Nx_base`**: Nodos en dirección horizontal (aumentar → mayor resolución en X)
- **`Ny_base`**: Nodos en dirección vertical (aumentar → mayor resolución en Y)
- **`h_base`**: Paso espacial base (reducir → mayor precisión)
- **`blocks_base`**: Coordenadas de los bloques/obstáculos dentro del dominio

### 3. Configurar las Simulaciones

Cada simulación se crea con una instancia de `NewtonRaphson`. Localiza las instancias en `main.py` (líneas 135-175) y modifica:

```python
NewtonRaphson(
    Nx=nxA,              # Nodos en X (de get_scaled_config)
    Ny=nyA,              # Nodos en Y (de get_scaled_config)
    blocks=bloqA,        # Bloques (de get_scaled_config)
    V0=1,                # Velocidad característica para visualización
    v_ij=0.01,           # Valor inicial de velocidad en el dominio
    h=h_valA,            # Paso espacial (de get_scaled_config)
    v0ParedSuperior=1    # Velocidad en la pared superior (0 o 1)
)
```

**Parámetros clave:**

- **`V0`**: Escala de velocidad característica (típicamente 1 o 2)
- **`v_ij`**: Valor inicial para todos los nodos internos (ej: 0.01, 0.05)
- **`v0ParedSuperior`**: Condición de frontera en la pared superior (0 = sin movimiento, 1 = con movimiento)

### 4. Habilitar/Deshabilitar Simulaciones

Al final de `main.py` (líneas 310-315), controla qué simulaciones ejecutar:

```python
generateSimulationA(True)   # True = ejecutar, False = omitir
generateSimulationB(False, "JacobianoInversa")  # False = omitir (ahorra tiempo)
```

**Nota:** La simulación B con `FACTORB = 3` puede tardar varios minutos. Ponla en `False` para pruebas rápidas.

### 5. Ajustar Criterios de Convergencia

En la función `iteration()` (líneas 51-54), puedes modificar:

```python
max_iterations = 200    # Máximo de iteraciones permitidas
epsilon = 1e-13         # Tolerancia para el residuo
delta = 1e-13           # Tolerancia para el cambio entre iteraciones
```

### 6. Activar Pruebas de Propiedades de la Matriz en Cada Iteración

Para ver las propiedades de la matriz Jacobiana durante **cada iteración** del método de Newton-Raphson, debes pasar el parámetro `proves=True` al llamar la función `iteration()`.

**Qué información muestra `proves=True`:**

- Número de condición del Jacobiano
- Si cumple con la condición de Richardson
- Si la matriz es simétrica
- Si es diagonalmente dominante
- Si es de rango completo

**Ejemplo de uso:**

```python
# Activar pruebas en una simulación específica
iteration(NewtonRaphsonItGCA, "gradiente Conjugado", saveHistory=True, proves=True)
```

**Para activar en todas las simulaciones del proyecto:**

Localiza las llamadas a `iteration()` en `main.py` y agrega `proves=True`:

```python
# Línea ~188 - Comparación de métodos
iteration(NewtonRaphsonItGCA, "gradiente Conjugado", proves=True)
iteration(NewtonRapshonItJIA, "jacobiano Inverso", proves=True)

# Dentro de generateSimulationA()
iteration(NewtonRaphsonItGCA2, "gradiente Conjugado", proves=True)
iteration(NewtonRaphsonItGCA3, "gradiente Conjugado", proves=True)
historial_pasosA = iteration(NewtonRaphsonItGCA, "gradiente Conjugado", True, proves=True)

# Dentro de generateSimulationB()
iteration(NewtonRaphsonItGCB, metodo, proves=True)
iteration(NewtonRaphsonItGCB2, metodo, proves=True)
iteration(NewtonRaphsonItGCB3, metodo, proves=True)
historial_pasosB = iteration(NewtonRaphsonItGCB, metodo, True, proves=True)
```

**Nota:** Activar `proves=True` aumentará el tiempo de ejecución, especialmente en mallas grandes (factor > 2), ya que calcula propiedades adicionales de la matriz en cada iteración.

---

## Detalles Técnicos de `get_scaled_config(factor)`

**Firma:** `get_scaled_config(factor: float) -> (Nx_new, Ny_new, h_new, blocks_new)`

**Funcionamiento:**

- Parte de una configuración base: `Nx_base = 81`, `Ny_base = 9`, `h_base = 5.0`, `blocks_base = [[10, 20, 0, 2], [58, 80, 5, 8]]`
- Calcula `Nx_new = int(Nx_base * factor)` y `Ny_new = int(Ny_base * factor)`
- Ajusta el paso espacial a `h_new = h_base / factor` (factor > 1 → mayor resolución)
- Escala las coordenadas de los bloques y corrige límites
- Devuelve `(Nx_new, Ny_new, h_new, blocks_new)`

**Ejemplos:**

- `get_scaled_config(1)` → Malla base 81×9
- `get_scaled_config(3)` → Malla ampliada 243×27
- `get_scaled_config(0.5)` → Malla reducida 40×4

---

## Ejemplos de Uso Común

### Ejemplo 1: Prueba Rápida

```python
# En main.py, líneas 125-126
FACTORA = 0.5  # Malla pequeña para pruebas rápidas
FACTORB = 1    # Reducir para pruebas más rápidas

# Al final del archivo
generateSimulationA(True)
generateSimulationB(False, "gradiente Conjugado")  # Desactivar B
```

### Ejemplo 2: Simulación de Alta Precisión

```python
# En main.py
FACTORA = 2
FACTORB = 5  # Advertencia: puede tardar mucho tiempo

# Dentro de get_scaled_config
h_base = 2.5  # Paso espacial más fino
```

### Ejemplo 3: Comparar Métodos de Resolución

```python
# Ejecutar con gradiente conjugado
iteration(NewtonRaphsonItGCA, "gradiente Conjugado", proves=True)

# Ejecutar con jacobiano inverso
iteration(NewtonRapshonItJIA, "jacobiano Inverso", proves=True)
```

---

## Solución de Problemas

- **El programa tarda mucho:** Reduce `FACTORB` o ponlo en `False` en `generateSimulationB(False, ...)`
- **Errores de memoria:** Reduce `Nx_base` y `Ny_base` en `get_scaled_config`
- **No converge:** Ajusta `epsilon`, `delta` o aumenta `max_iterations` en `iteration()`
- **Quiero ver las propiedades de la matriz:** Pasa `proves=True` a la función `iteration()`

---

## Estructura del Flujo de Trabajo

1. **Configuración** → Modifica `FACTORA`, `FACTORB` y parámetros en `get_scaled_config`
2. **Ejecución** → Corre `python main.py`
3. **Visualización** → Navega por las diapositivas en `SimulationViewer` usando los botones
4. **Análisis** → Observa convergencia, mapas de calor, splines y animaciones
