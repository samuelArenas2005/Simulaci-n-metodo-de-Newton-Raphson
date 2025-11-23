import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RectBivariateSpline

class SimulationViewer:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Visor de Simulación CFD - Newton Raphson")
        self.root.geometry("1100x800")

        # Estado interno
        self.slides = []
        self.current_index = 0
        self.current_anim = None # Referencia vital para que no se congele la animación

        # --- Interfaz Gráfica ---
        
        # 1. Título y Contador
        self.header_frame = tk.Frame(self.root)
        self.header_frame.pack(side=tk.TOP, fill=tk.X, pady=5)
        
        self.lbl_title = tk.Label(self.header_frame, text="", font=("Arial", 14, "bold"))
        self.lbl_title.pack()
        
        self.lbl_counter = tk.Label(self.header_frame, text="0 / 0", font=("Arial", 10), fg="gray")
        self.lbl_counter.pack()

        # 2. Área del Gráfico (Matplotlib embed)
        # Creamos una figura vacía que pasaremos a tus funciones para que pinten en ella
        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Barra de herramientas (Zoom, Pan, Guardar)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.root)
        self.toolbar.update()

        # 3. Controles de Navegación
        btn_frame = tk.Frame(self.root, pady=10, bg="#f0f0f0")
        btn_frame.pack(side=tk.BOTTOM, fill=tk.X)

        tk.Button(btn_frame, text="<< Anterior", command=self.prev_slide, width=15).pack(side=tk.LEFT, padx=20)
        tk.Button(btn_frame, text="Siguiente >>", command=self.next_slide, width=15).pack(side=tk.RIGHT, padx=20)

    def add_slide(self, title, draw_function, *args):
        """
        Registra una vista.
        draw_function: Función que recibe (fig, ax, *args) y dibuja.
        """
        self.slides.append({
            'title': title,
            'func': draw_function,
            'args': args
        })

    def update_view(self):
        if not self.slides: return

        # 1. Limpieza
        self.ax.clear()
        # Si hay una animación corriendo, la detenemos
        if self.current_anim:
            self.current_anim.event_source.stop()
            self.current_anim = None

        # 2. Obtener datos actuales
        slide = self.slides[self.current_index]
        title = slide['title']
        func = slide['func']
        args = slide['args']

        # 3. Actualizar textos
        self.lbl_title.config(text=title)
        self.lbl_counter.config(text=f"Vista {self.current_index + 1} de {len(self.slides)}")

        # 4. EJECUTAR LA FUNCIÓN DE DIBUJO
        # Aquí ocurre la magia: llamamos a la función adaptada pasando NUESTRA figura y eje
        result = func(self.fig, self.ax, *args)

        # 5. Si devolvió una animación, guardarla
        if isinstance(result, animation.FuncAnimation):
            self.current_anim = result

        # 6. Refrescar lienzo
        self.canvas.draw()

    def next_slide(self):
        if self.current_index < len(self.slides) - 1:
            self.current_index += 1
            self.update_view()

    def prev_slide(self):
        if self.current_index > 0:
            self.current_index -= 1
            self.update_view()

    def run(self):
        self.update_view()
        self.root.mainloop()


def view_heatmap(fig, ax, plot_obj):
    """Adaptador para showPlot"""
    NX, NY = plot_obj.Nx, plot_obj.Ny
    matriz = np.array(plot_obj.vec).reshape(NY, NX)
    
    im = ax.imshow(matriz, cmap=plt.cm.viridis, vmin=0.0001, vmax=plot_obj.V0, origin='lower')
    ax.set_title(f"Mapa de Calor Original")
    ax.set_xlabel("i (x)")
    ax.set_ylabel("j (y)")
    # Nota: En animaciones complejas las colorbars repetidas pueden dar error, 
    # aquí simplificamos asumiendo que el toolbar permite ver valores.
    return im

def view_spline_static(fig, ax, plot_obj, suavizado):
    """Adaptador para showPlotSpline"""
    NX, NY = plot_obj.Nx, plot_obj.Ny
    
    # --- Tu Lógica de Interpolación ---
    x = np.arange(NX)
    y = np.arange(NY)
    z = np.array(plot_obj.vec).reshape(NY, NX)
    spline = RectBivariateSpline(y, x, z)

    x_new = np.linspace(0, NX - 1, NX * suavizado)
    y_new = np.linspace(0, NY - 1, NY * suavizado)
    z_new = spline(y_new, x_new)
    
    # --- Tu Lógica de Máscara ---
    X_grid, Y_grid = np.meshgrid(x_new, y_new)
    if hasattr(plot_obj, 'blocks'):
        for block in plot_obj.blocks:
            x1, x2, y1, y2 = block
            mask = (X_grid >= x1) & (X_grid <= x2) & (Y_grid >= y1) & (Y_grid <= y2)
            z_new[mask] = 0.0
    z_new[z_new < 0] = 0 

    im = ax.imshow(z_new, cmap=plt.cm.viridis, vmin=0.0001, vmax=plot_obj.V0, 
                   origin='lower', extent=[0, NX - 1, 0, NY - 1])
    ax.set_title(f"Spline Suavizado x{suavizado}")
    return im

def view_anim_history(fig, ax, plot_obj, historial, intervalo):
    """Adaptador para animateHistoryLoop"""
    NX, NY = plot_obj.Nx, plot_obj.Ny
    
    matriz_inicial = historial[0].reshape(NY, NX)
    im = ax.imshow(matriz_inicial, cmap=plt.cm.viridis, vmin=0.0001, vmax=plot_obj.V0, origin='lower')
    ax.set_xlabel("i (x)")
    ax.set_ylabel("j (y)")
    
    def update(frame):
        vec_actual = historial[frame]
        im.set_data(vec_actual.reshape(NY, NX))
        ax.set_title(f"Iteración: {frame}")
        return [im]

    # Importante: Usamos 'fig' que viene del visor
    ani = animation.FuncAnimation(fig, update, frames=len(historial), 
                                  interval=intervalo, blit=False, repeat=True)
    return ani

def view_anim_spline(fig, ax, plot_obj, max_suavizado, intervalo):
    """Adaptador para animateSpline"""
    NX, NY = plot_obj.Nx, plot_obj.Ny
    x_orig = np.arange(NX)
    y_orig = np.arange(NY)
    z_orig = np.array(plot_obj.vec).reshape(NY, NX)
    spline = RectBivariateSpline(y_orig, x_orig, z_orig)
    
    # Configuración estática inicial
    im = ax.imshow(np.zeros((NY, NX)), cmap=plt.cm.viridis, vmin=0.0001, vmax=plot_obj.V0, 
                   origin='lower', extent=[0, NX - 1, 0, NY - 1])

    def update(frame):
        factor = frame + 1
        # --- Tu lógica de cálculo por cuadro ---
        x_new = np.linspace(0, NX - 1, NX * factor)
        y_new = np.linspace(0, NY - 1, NY * factor)
        z_new = spline(y_new, x_new)
        
        X_grid, Y_grid = np.meshgrid(x_new, y_new)
        if hasattr(plot_obj, 'blocks'):
            for block in plot_obj.blocks:
                x1, x2, y1, y2 = block
                mask = (X_grid >= x1) & (X_grid <= x2) & (Y_grid >= y1) & (Y_grid <= y2)
                z_new[mask] = 0.0
        z_new[z_new < 0] = 0 
        
        # OJO: Para cambiar la matriz de datos en imshow cuando cambia el tamaño,
        # lo mejor es limpiar y redibujar o resetear data y extent.
        # Aquí usamos set_data + set_extent
        im.set_data(z_new)
        im.set_extent([0, NX - 1, 0, NY - 1])
        ax.set_title(f"Spline Animado: Resolución x{factor}")
        return [im]

    ani = animation.FuncAnimation(fig, update, frames=max_suavizado, 
                                  interval=intervalo, blit=False)
    return ani