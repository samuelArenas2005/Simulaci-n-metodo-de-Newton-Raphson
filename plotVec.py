import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.animation as animation

class Plot:

    #Valores iniciales de nuestra simulación que permite 
    def __init__(self,vector):
        self.vec = vector.vec
        self.Nx = vector.Nx
        self.Ny = vector.Ny
        self.V0 = vector.V0
        self.blocks = vector.blocks
    
    #Muestra el mapa de calor del vector de velocidades
    def showPlot(self, condicion=True):
        if condicion:
            x0 = self.vec
            # --- visualización ---
            # Usar las dimensiones de la clase en lugar de valores hardcodeados
            NX, NY = self.Nx, self.Ny
            matriz = np.array(x0).reshape(NY, NX)

            # crear un mapa de calor con tamaño proporcional a las dimensiones
            fig_width = max(10, NX * 0.15)  # Ancho mínimo 10, escalado con NX
            fig_height = max(3, NY * 0.5)   # Alto mínimo 3, escalado con NY
            plt.figure(figsize=(fig_width, fig_height))
            
            cmap = plt.cm.viridis  # colormap principal

            # plot con imshow
            im = plt.imshow(matriz[::-1], cmap=cmap, vmin=0.0001, vmax=self.V0)  
            # matriz[::-1] → para que se vea de abajo hacia arriba (como ejes cartesianos)
            # vmax usa V0 en lugar de 1 hardcodeado

            plt.colorbar(im, label="u valores")
            plt.title(f"Distribución de valores en la malla ({NX}x{NY})")
            plt.xlabel("i (x)")
            plt.ylabel("j (y)")
            plt.show()
    
    #Metodo que muestra la grafica en cada iteración
    def animateHistoryLoop(self, historial, intervalo=200):
        """
        Reproduce el historial de iteraciones en bucle infinito.
        
        historial: Lista de vectores obtenida de iteration_con_historial
        intervalo: Tiempo en ms entre cuadros (ej. 500 = medio segundo, lento).
        """
        NX, NY = self.Nx, self.Ny
        total_frames = len(historial)
        
        # Configurar figura
        fig_width = max(10, NX * 0.15)
        fig_height = max(3, NY * 0.5)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        # Configurar imagen inicial
        # Tomamos el primer cuadro para configurar la escala de colores
        matriz_inicial = historial[0].reshape(NY, NX)
        
        cmap = plt.cm.viridis
        im = ax.imshow(matriz_inicial, 
                       cmap=cmap, 
                       vmin=0.0001, 
                       vmax=self.V0, 
                       origin='lower') # Importante: origin lower
        
        plt.colorbar(im, label="u valores")
        ax.set_xlabel("i (x)")
        ax.set_ylabel("j (y)")
        title = ax.set_title("Iniciando...")

        # Función de actualización para la animación
        def update(frame):
            # 1. Obtener datos del historial
            vec_actual = historial[frame]
            matriz = vec_actual.reshape(NY, NX)
            
            # 2. Actualizar imagen
            im.set_data(matriz)
            
            # 3. Actualizar título
            title.set_text(f"Iteración: {frame}/{total_frames-1}")
            
            return [im, title]

        # Crear la animación
        # repeat=True hace que sea un BUCLE
        # repeat_delay=1000 espera 1 segundo al final antes de volver a empezar
        ani = animation.FuncAnimation(fig, update, frames=total_frames, 
                                      interval=intervalo, blit=False, 
                                      repeat=True, repeat_delay=1000)
        
        print("Reproduciendo animación en bucle (Cierre la ventana para salir)...")
        plt.show()

    #Muestra el mapa de calor del vector de velocidades con sus valores exactos
    def showPlotDetail(self, condicion=False):
      if condicion:
          x0 = self.vec
          NX, NY = self.Nx, self.Ny
          matriz = np.array(x0).reshape(NY, NX)

          escala = 0.4
          plt.figure(figsize=(NX * escala, NY * escala))

          cmap = plt.cm.viridis
          im = plt.imshow(matriz, cmap=cmap, vmin=0.0001, vmax=1, origin="lower")

          plt.colorbar(im, label="u valores")
          plt.title("Distribución de valores finales en la malla")
          plt.xlabel("i (x)")
          plt.ylabel("j (y)")

          # cuadrícula
          ax = plt.gca()
          ax.set_xticks(np.arange(-0.5, NX, 1), minor=True)
          ax.set_yticks(np.arange(-0.5, NY, 1), minor=True)
          ax.grid(which="minor", color="white", linestyle='-', linewidth=0.5)
          ax.tick_params(which="minor", bottom=False, left=False)

          fontsize = 8
          for i in range(NY):
              for j in range(NX):
                  valor = matriz[i, j]

                  val_trunc = np.trunc(valor * 100.0) / 100.0

                  # si es cero tras truncar -> no mostrar nada
                  if abs(val_trunc) < 1e-12:
                      continue

                  # extraer la parte "xx" correspondiente a los centésimos:
                  # por ejemplo 1.237 -> 237 -> 237 % 100 = 37 -> ".37"
                  centesimos = int(abs(np.trunc(val_trunc * 100.0))) % 100

                  # si la parte fraccional es 0 -> no mostrar nada (eran enteros .00)
                  if centesimos == 0:
                      continue

                  sign = "-" if val_trunc < 0 else ""
                  texto = f"{sign}.{centesimos:02d}"

                  # color del texto (blanco sobre el mapa)
                  color_texto = "white"

                  plt.text(j, i, texto, ha="center", va="center",
                          color=color_texto, fontsize=fontsize)

          plt.show()

    #Muestra la malla de calor con los indices que le corresponde, lo cual nos permite guiarnos más fácil para los calculos
    def showPlotIndices(self, condicion=False):
      if condicion:
          x0 = self.vec
          NX, NY = self.Nx, self.Ny
          matriz = np.array(x0).reshape(NY, NX)

          escala = 0.4
          plt.figure(figsize=(NX * escala, NY * escala)) # figura grande para que se lea bien
          cmap = plt.cm.viridis

          im = plt.imshow(matriz, cmap=cmap, vmin=0.0001, vmax=1, origin="lower")

          plt.colorbar(im, label="u valores")
          plt.title("Distribución de valores en la malla con índices")
          plt.xlabel("i (x)")
          plt.ylabel("j (y)")

          # Dibujar bordes de cada celda
          for i in range(NX + 1):
              plt.axvline(i - 0.5, color='gray', linewidth=0.3)
          for j in range(NY + 1):
              plt.axhline(j - 0.5, color='gray', linewidth=0.3)

          # Colocar el índice en cada celda
          for j in range(NY):
              for i in range(NX):
                  idx = j * NX + i
                  plt.text(i, j, str(idx),
                          ha='center', va='center',
                          color='white', fontsize=6)

          # ticks que correspondan exactamente a los cuadros
          plt.xticks(np.arange(NX), np.arange(NX))
          plt.yticks(np.arange(NY), np.arange(NY))
          
          plt.show()
          


    # Muestra el mapa de calor suavizado usando interpolación Spline Bicúbica
    def showPlotSpline(self, suavizado=10):
        NX, NY = self.Nx, self.Ny
        
        # 1. Interpolación (Igual que antes)
        x = np.arange(NX)
        y = np.arange(NY)
        z = np.array(self.vec).reshape(NY, NX)

        spline = RectBivariateSpline(y, x, z)

        x_new = np.linspace(0, NX - 1, NX * suavizado)
        y_new = np.linspace(0, NY - 1, NY * suavizado)
        z_new = spline(y_new, x_new)
        
        # Creamos matrices de coordenadas X e Y para cada píxel de la nueva imagen suave
        # Esto nos dice en qué coordenada "real" está cada píxel interpolado
        X_grid, Y_grid = np.meshgrid(x_new, y_new)

        # Recorremos cada bloque definido en tu simulación
        if hasattr(self, 'blocks'):
            for block in self.blocks:
                x1, x2, y1, y2 = block
                
                # Creamos una máscara booleana: ¿Qué píxeles caen dentro de este bloque?
                # Usamos los límites del bloque. 
                # Nota: Dependiendo de tu lógica de índices, quizás necesites ajustar +/- 0.5
                mask = (X_grid >= x1) & (X_grid <= x2) & (Y_grid >= y1) & (Y_grid <= y2)
                
                # Forzamos esos píxeles a CERO absoluto
                z_new[mask] = 0.0

        # Limpieza adicional para oscilaciones negativas fuera de los bloques
        z_new[z_new < 0] = 0 

        # ---------------------------------------------------------
        # 3. Visualización
        # ---------------------------------------------------------
        fig_width = max(10, NX * 0.15)
        fig_height = max(3, NY * 0.5)
        plt.figure(figsize=(fig_width, fig_height))
        
        cmap = plt.cm.viridis
        
        im = plt.imshow(z_new, 
                        cmap=cmap, 
                        vmin=0.0001, 
                        vmax=self.V0, 
                        origin='lower',
                        extent=[0, NX - 1, 0, NY - 1]) 

        plt.colorbar(im, label="u valores (Interpolado)")
        plt.title(f"Distribución Suavizada con Obstáculos ({NX}x{NY})")
        plt.xlabel("i (x)")
        plt.ylabel("j (y)")
        plt.show()
    
    def animateSpline(self, max_suavizado=10, intervalo=1000):
        """
        Animación de la interpolación.
        
        intervalo: Tiempo en milisegundos entre cuadros (1000ms = 1 segundo).
                   Cuanto mayor sea este número, más lenta será la animación.
        """
        NX, NY = self.Nx, self.Ny
        
        # 1. Pre-calcular (Igual que antes)
        x_orig = np.arange(NX)
        y_orig = np.arange(NY)
        z_orig = np.array(self.vec).reshape(NY, NX)
        
        spline = RectBivariateSpline(y_orig, x_orig, z_orig)

        # 2. Configurar figura
        fig_width = max(10, NX * 0.15)
        fig_height = max(3, NY * 0.5)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        cmap = plt.cm.viridis

        def update(frame):
            factor = frame + 1
            
            # A. Generar nueva malla
            x_new = np.linspace(0, NX - 1, NX * factor)
            y_new = np.linspace(0, NY - 1, NY * factor)
            
            # B. Calcular valores interpolados
            z_new = spline(y_new, x_new)
            
            # C. Enmascarar Bloques (Ponerlos a cero estricto)
            X_grid, Y_grid = np.meshgrid(x_new, y_new)
            
            if hasattr(self, 'blocks'):
                for block in self.blocks:
                    x1, x2, y1, y2 = block
                    mask = (X_grid >= x1) & (X_grid <= x2) & (Y_grid >= y1) & (Y_grid <= y2)
                    z_new[mask] = 0.0
            
            z_new[z_new < 0] = 0 

            # D. Graficar
            ax.clear()
            
            im = ax.imshow(z_new, 
                           cmap=cmap, 
                           vmin=0.0001, 
                           vmax=self.V0, 
                           origin='lower',
                           extent=[0, NX - 1, 0, NY - 1],
                           interpolation='nearest')
            
            ax.set_title(f"Resolución x{factor} (Interpolación Spline)")
            ax.set_xlabel("i (x)")
            ax.set_ylabel("j (y)")
            
            # --- NOTA: Se eliminó el código que dibujaba los bordes blancos ---

            return [im]

        # 3. Ejecutar animación
        # intervalo controla la velocidad (en milisegundos)
        ani = animation.FuncAnimation(fig, update, frames=max_suavizado, interval=intervalo, blit=False)
        
        # Barra de color estática
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0.0001, vmax=self.V0))
        fig.colorbar(sm, ax=ax, label="u valores")

        plt.show()