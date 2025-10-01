# sim_gui.py
import tkinter as tk
from tkinter import ttk, messagebox
from tkinter.scrolledtext import ScrolledText

import numpy as np
import random
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# ---------------------------
# Clase Vector (robusta)
# ---------------------------
class Vector:
    def __init__(self, NX=81, NY=9, vij=0.001, blockA=None, blockB=None, cut_x_physical=None, Lx=400):
        self.NX = int(NX)
        self.NY = int(NY)
        self.n = self.NX * self.NY
        self.vec = np.zeros(self.n)
        self.vectFunction = np.zeros(self.n)
        self.matrixJacobiana = np.zeros((self.n, self.n))
        self.vij = float(vij)
        self.blockA = blockA
        self.blockB = blockB
        self.Lx = float(Lx)
        # cut_i proporcional a posición física de 58 (mantener compatibilidad)
        if cut_x_physical is None:
            self.cut_i = int((58.0 / 400.0) * max(1, self.NX))
        else:
            # si te dan cut_x_physical en unidades físicas, lo convierto
            self.cut_i = int((cut_x_physical / self.Lx) * self.NX)

    def _in_block(self, i, j, block):
        if block is None:
            return False
        i0,i1,j0,j1 = block
        return (i >= i0 and i <= i1 and j >= j0 and j <= j1)

    def vecInicial(self, velocidad_init=1, post_block_recover=False):
        V0 = velocidad_init
        rango = [1,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0]
        x0 = self.vec
        k = 0
        cols_per_band = max(1, int(self.NX / 10))

        # blockB horizontal extents (si existe)
        if self.blockB:
            b_i0, b_i1, b_j0, b_j1 = self.blockB
        else:
            b_i0 = None
            b_i1 = None

        for i in range(self.NX):
            if (i % cols_per_band == cols_per_band - 1) and (k + 1 < len(rango)):
                k = min(k + 1, len(rango)-1)
            for j in range(self.NY):
                idx = j*self.NX + i
                # izq vertical
                if (i == 0 and j > 0):
                    x0[idx] = V0
                    continue

                # borde superior (j == NY-1) dinámica con respecto a blockB
                if j == (self.NY - 1):
                    if b_i0 is not None:
                        if i < b_i0:
                            x0[idx] = V0
                        elif b_i0 <= i <= b_i1:
                            x0[idx] = 0.0
                        else:
                            # comportamiento post-block: por defecto 0, pero si post_block_recover True -> V0
                            x0[idx] = V0 if post_block_recover else 0.0
                        continue
                    else:
                        # sin blockB definido: usar cut_i (comportamiento legado)
                        if i < self.cut_i:
                            x0[idx] = V0
                        else:
                            x0[idx] = 0.0
                        continue

                # bordes y bloques
                if (j == 0 or i == self.NX - 1 or self._in_block(i, j, self.blockB) or self._in_block(i, j, self.blockA)):
                    x0[idx] = 0.0
                    continue

                # caso especial similar al original (siempre tratar bloque A original)
                if (not self.blockA) and (10 <= i <= 20) and (1 <= j <= 2):
                    x0[idx] = random.uniform(rango[9], rango[10])
                    continue

                x0[idx] = random.uniform(rango[k], rango[min(k+1, len(rango)-1)])

    def cal_function(self, post_block_recover=False):
        val_Vij = self.vij
        x0 = self.vec
        NX = self.NX
        NY = self.NY
        n = NX * NY

        if self.vectFunction.size != n:
            self.vectFunction = np.zeros(n)
        self.vectFunction.fill(0.0)

        # blockB horizontal extents (si existe)
        if self.blockB:
            b_i0, b_i1, b_j0, b_j1 = self.blockB
        else:
            b_i0 = None
            b_i1 = None

        for j in range(NY):
            for i in range(NX):
                idx = j*NX + i
                if idx < 0 or idx >= n:
                    continue

                # izquierda vertical
                if (i == 0 and j > 0):
                    self.vectFunction[idx] = 0.0
                    continue

                # borde superior dinámico
                if j == (NY - 1):
                    if b_i0 is not None:
                        # todas las condiciones del borde se tratan como 0 en vectFunction
                        self.vectFunction[idx] = 0.0
                        continue
                    else:
                        self.vectFunction[idx] = 0.0
                        continue

                # otras condiciones de borde / bloques
                if (j == 0 or i == NX - 1 or self._in_block(i, j, self.blockB) or self._in_block(i, j, self.blockA)):
                    self.vectFunction[idx] = 0.0
                    continue

                # comprobar vecinos
                if not ( (i+1) < NX and (i-1) >= 0 and (j+1) < NY and (j-1) >= 0 ):
                    self.vectFunction[idx] = 0.0
                    continue

                ip = idx + 1
                im = idx - 1
                jp = idx + NX
                jm = idx - NX

                self.vectFunction[idx] = (1.0/4.0) * (
                    x0[ip] + x0[im] + x0[jp] + x0[jm]
                    - (5.0/2.0) * x0[idx] * (x0[ip] - x0[im])
                    - (5.0/2.0) * (val_Vij) * (x0[jp] - x0[jm])
                ) - x0[idx]

    def cal_jacobiano(self):
        v_ij = self.vij
        vec_0 = self.vec
        NX = self.NX
        NY = self.NY
        n = self.n

        self.matrixJacobiana = np.zeros((n, n))

        valores_derechos = [ecuacion * NX - 1 for ecuacion in range(1, NY+1)]
        valores_bloque_A = []
        valores_bloque_B = []
        if self.blockA:
            i0,i1,j0,j1 = self.blockA
            valores_bloque_A = [k * NX + j for j in range(i0, i1+1) for k in range(j0, j1+1)]
        if self.blockB:
            i0,i1,j0,j1 = self.blockB
            valores_bloque_B = [k * NX + j for j in range(i0, i1+1) for k in range(j0, j1+1)]

        for ecuacion in range(n):
            if ((ecuacion % NX == 0) or (ecuacion < NX) or (ecuacion >= NX*(NY-1))
                    or ecuacion in valores_derechos or ecuacion in valores_bloque_A or ecuacion in valores_bloque_B):
                self.matrixJacobiana[ecuacion, ecuacion] = 1.0
                continue

            i = ecuacion
            if (i+1) >= n or (i-1) < 0 or (i+NX) >= n or (i-NX) < 0:
                self.matrixJacobiana[ecuacion, ecuacion] = 1.0
                continue

            try:
                self.matrixJacobiana[i, i]   = -1.0 - (5.0/8.0) * (vec_0[i+1] - vec_0[i-1])
                self.matrixJacobiana[i, i+1] =  1.0/4.0 - (5.0/8.0) * vec_0[i]
                self.matrixJacobiana[i, i-1] =  1.0/4.0 + (5.0/8.0) * vec_0[i]
                self.matrixJacobiana[i, i+NX] = 1.0/4.0 - (5.0/8.0) * v_ij
                self.matrixJacobiana[i, i-NX] = 1.0/4.0 + (5.0/8.0) * v_ij
            except Exception:
                self.matrixJacobiana[ecuacion, ecuacion] = 1.0

    def newVector(self):
        try:
            delta = np.linalg.solve(self.matrixJacobiana, self.vectFunction)
        except Exception:
            delta = np.linalg.pinv(self.matrixJacobiana).dot(self.vectFunction)
        self.vec = self.vec - delta

    # figuras: devuelven objetos matplotlib.Figure
    def fig_showPlotDetail(self):
        x0 = self.vec
        NX, NY = self.NX, self.NY
        matriz = np.array(x0).reshape(NY, NX)
        escala = 0.25
        fig = plt.Figure(figsize=(12, 6))
        ax = fig.add_subplot(111)
        im = ax.imshow(matriz, cmap=plt.cm.viridis, vmin=0.0001, vmax=1, origin="lower")
        fig.colorbar(im, ax=ax, label="u valores")
        ax.set_title("Distribución valores (detalle)")
        ax.set_xlabel("i (x)")
        ax.set_ylabel("j (y)")
        ax.set_xticks(np.arange(0, NX, max(1, NX//10)))
        ax.set_yticks(np.arange(0, NY, max(1, NY//4)))
        # cuadricula y valores
        ax.set_xticks(np.arange(-0.5, NX, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, NY, 1), minor=True)
        """   ax.grid(which="minor", color="white", linestyle='-', linewidth=0.5) """
        ax.tick_params(which="minor", bottom=False, left=False)
        
        # marcar bloque B (si existe)
       
        return fig

# ---------------------------
# GUI
# ---------------------------
class App:
    def __init__(self, root):
        self.root = root
        root.title("Simulación - Interfaz completa")
        frm = ttk.Frame(root, padding=8)
        frm.grid(row=0, column=0, sticky="nsew")
        row = 0

        # Física / mallado
        ttk.Label(frm, text="Lx (long. física)").grid(row=row, column=0, sticky="w")
        self.e_Lx = ttk.Entry(frm); self.e_Lx.insert(0,"400"); self.e_Lx.grid(row=row, column=1, sticky="ew")
        ttk.Label(frm, text="Ly (long. física)").grid(row=row, column=2, sticky="w")
        self.e_Ly = ttk.Entry(frm); self.e_Ly.insert(0,"40"); self.e_Ly.grid(row=row, column=3, sticky="ew")
        row += 1

        ttk.Label(frm, text="hx (paso x)").grid(row=row, column=0, sticky="w")
        self.e_hx = ttk.Entry(frm); self.e_hx.insert(0,"5"); self.e_hx.grid(row=row, column=1, sticky="ew")
        ttk.Label(frm, text="hy (paso y)").grid(row=row, column=2, sticky="w")
        self.e_hy = ttk.Entry(frm); self.e_hy.insert(0,"5"); self.e_hy.grid(row=row, column=3, sticky="ew")
        row += 1

        ttk.Label(frm, text="NX (nodos x)").grid(row=row, column=0, sticky="w")
        self.var_nx = tk.StringVar(value="81")
        self.l_nx = ttk.Entry(frm, textvariable=self.var_nx, state="readonly"); self.l_nx.grid(row=row, column=1, sticky="ew")
        ttk.Label(frm, text="NY (nodos y)").grid(row=row, column=2, sticky="w")
        self.var_ny = tk.StringVar(value="9")
        self.l_ny = ttk.Entry(frm, textvariable=self.var_ny, state="readonly"); self.l_ny.grid(row=row, column=3, sticky="ew")
        row += 1

        # parámetros iniciales
        ttk.Label(frm, text="V0 (vel init)").grid(row=row, column=0, sticky="w")
        self.e_v0 = ttk.Entry(frm); self.e_v0.insert(0,"1"); self.e_v0.grid(row=row, column=1, sticky="ew")
        ttk.Label(frm, text="v_ij").grid(row=row, column=2, sticky="w")
        self.e_vij = ttk.Entry(frm); self.e_vij.insert(0,"0.001"); self.e_vij.grid(row=row, column=3, sticky="ew")
        row += 1

        # bloques
        ttk.Label(frm, text="Block A (i0-i1,j0-j1 ó x0-x1,y0-y1)").grid(row=row, column=0, sticky="w")
        self.e_blockA = ttk.Entry(frm); self.e_blockA.insert(0,"10-20,1-2"); self.e_blockA.grid(row=row, column=1, sticky="ew")
        ttk.Label(frm, text="Block B (i0-i1,j0-j1 ó x0-x1,y0-y1)").grid(row=row, column=2, sticky="w")
        self.e_blockB = ttk.Entry(frm); self.e_blockB.insert(0,"58-80,5-8"); self.e_blockB.grid(row=row, column=3, sticky="ew")
        row += 1

        # opciones y convergencia
        self.post_block_recover_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(frm, text="Después del bloque superior: volver a V0", variable=self.post_block_recover_var).grid(row=row, column=0, columnspan=2, sticky="w")
        self.show_mini_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(frm, text="Mostrar mini-plots (mapa + índices)", variable=self.show_mini_var).grid(row=row, column=2, columnspan=2, sticky="w")
        row += 1

        ttk.Label(frm, text="Max iter").grid(row=row, column=0, sticky="w")
        self.e_maxiter = ttk.Entry(frm); self.e_maxiter.insert(0,"100"); self.e_maxiter.grid(row=row, column=1, sticky="ew")
        ttk.Label(frm, text="Tolerancia").grid(row=row, column=2, sticky="w")
        self.e_tol = ttk.Entry(frm); self.e_tol.insert(0,"1e-10"); self.e_tol.grid(row=row, column=3, sticky="ew")
        row += 1

        # botones
        self.btn_compute = ttk.Button(frm, text="Calcular y Mostrar (detalle grande)", command=self.on_compute)
        self.btn_compute.grid(row=row, column=0, columnspan=2, sticky="ew", pady=6)
        self.btn_quit = ttk.Button(frm, text="Salir", command=root.quit)
        self.btn_quit.grid(row=row, column=2, columnspan=2, sticky="ew")
        row += 1

        # consola
        ttk.Label(frm, text="Consola").grid(row=row, column=0, sticky="w")
        self.console = ScrolledText(frm, height=10)
        self.console.grid(row=row+1, column=0, columnspan=4, sticky="nsew", pady=6)
        frm.columnconfigure(1, weight=1)
        frm.rowconfigure(row+1, weight=1)

    def log(self, *args):
        texto = " ".join(str(a) for a in args) + "\n"
        self.console.insert(tk.END, texto)
        self.console.see(tk.END)

    def parse_block(self, text, hx, hy, Lx, Ly):
        """
        Acepta 'i0-i1,j0-j1' (índices) o 'x0-x1,y0-y1' (coordenadas físicas).
        Devuelve (i0,i1,j0,j1) en índices de nodos o None si error.
        """
        try:
            parts = text.split(",")
            a = parts[0].strip()
            b = parts[1].strip()
            p0,p1 = [float(x) for x in a.split("-")]
            q0,q1 = [float(x) for x in b.split("-")]

            # detecto si son coordenadas físicas (si exceden dimensiones físicas o rango grande)
            is_physical = (p1 > Lx) or (q1 > Ly) or (p1 - p0 > Lx/2) or (q1 - q0 > Ly/2)

            if is_physical:
                i0 = max(0, int(p0 / hx))
                i1 = min(int(round(Lx / hx)), int(p1 / hx))
                j0 = max(0, int(q0 / hy))
                j1 = min(int(round(Ly / hy)), int(q1 / hy))
                return (i0, i1, j0, j1)
            else:
                i0 = int(p0); i1 = int(p1); j0 = int(q0); j1 = int(q1)
                return (i0, i1, j0, j1)
        except Exception:
            return None

    def on_compute(self):
        try:
            Lx = float(self.e_Lx.get())
            Ly = float(self.e_Ly.get())
            hx = float(self.e_hx.get())
            hy = float(self.e_hy.get())
            V0 = float(self.e_v0.get())
            vij = float(self.e_vij.get())
            max_iter = int(self.e_maxiter.get())
            tol = float(self.e_tol.get())
        except Exception as e:
            messagebox.showerror("Error", f"Parámetros numéricos inválidos: {e}")
            return

        # calcular NX, NY
        NX = int(round(Lx / hx)) + 1
        NY = int(round(Ly / hy)) + 1
        if NX < 3 or NY < 3:
            messagebox.showerror("Error", f"Nodos (NX,NY) demasiado pequeños: {NX},{NY}")
            return
        self.var_nx.set(str(NX)); self.var_ny.set(str(NY))

        # parsear bloques
        blockA = self.parse_block(self.e_blockA.get(), hx, hy, Lx, Ly)
        blockB = self.parse_block(self.e_blockB.get(), hx, hy, Lx, Ly)

        self.log("Parametros: Lx,Ly=", Lx, Ly, "hx,hy=", hx, hy, "NX,NY=", NX, NY)
        self.log("V0=", V0, "vij=", vij, "blockA=", blockA, "blockB=", blockB)
        self.log("max_iter=", max_iter, "tol=", tol, "post_block_recover=", self.post_block_recover_var.get())

        # crear Vector y ejecutar
        vec = Vector(NX=NX, NY=NY, vij=vij, blockA=blockA, blockB=blockB, Lx=Lx)
        vec.vecInicial(velocidad_init=V0, post_block_recover=self.post_block_recover_var.get())
        vec.cal_function(post_block_recover=self.post_block_recover_var.get())
        residuo_norm = np.linalg.norm(vec.vectFunction)
        self.log("Residuo inicial (norma):", residuo_norm)

        convergido = False
        for it in range(max_iter):
            self.log(" Iteración", it+1)
            vec.cal_function(post_block_recover=self.post_block_recover_var.get())
            residuo_norm = np.linalg.norm(vec.vectFunction)
            self.log("  Norma residuo:", residuo_norm)
            if residuo_norm < tol:
                self.log(f"¡Convergencia en {it+1} iteraciones!")
                convergido = True
                break
            vec.cal_jacobiano()
            vec.newVector()
        if not convergido:
            self.log("No convergió después de", max_iter, "iteraciones")

        # mostrar resultado (detalle grande por defecto)
        self.show_figures_window(vec)

    def show_figures_window(self, vec: Vector):
        win = tk.Toplevel(self.root)
        win.title("Resultados - Detalle grande")
        # Establecer tamaño inicial de ventana más grande
        win.geometry("1200x800")

        win.columnconfigure(0, weight=1); win.rowconfigure(0, weight=1)

        # Frame principal grande para el detalle
        f_main = ttk.Frame(win, padding=6)
        f_main.grid(row=0, column=0, sticky="nsew")

        # figura detalle grande - aumentar tamaño
        fig_detail = vec.fig_showPlotDetail()
        fig_detail.set_size_inches(16, 10)  # Tamaño más grande
        canvas_main = FigureCanvasTkAgg(fig_detail, master=f_main)
        canvas_main.draw()
        canvas_main.get_tk_widget().pack(fill="both", expand=True)

        # opcional: mini-plots en fila inferior
        if self.show_mini_var.get():
            win.rowconfigure(1, weight=0)
            f_row = ttk.Frame(win); f_row.grid(row=1, column=0, sticky="ew")
            f_map = ttk.Frame(f_row, padding=4); f_map.grid(row=0, column=0, sticky="nsew")
            f_idx = ttk.Frame(f_row, padding=4); f_idx.grid(row=0, column=1, sticky="nsew")
            fig_map = vec.fig_showPlot()      # mapa
            fig_idx = vec.fig_showPlotDetail()  # reuso detalle (o podrías llamar a showPlotIndices si la implementas)
            # reducir tamaño de miniaturas
            try:
                fig_map.set_size_inches(4,3)
                fig_idx.set_size_inches(4,3)
            except Exception:
                pass
            canvas_map = FigureCanvasTkAgg(fig_map, master=f_map); canvas_map.draw(); canvas_map.get_tk_widget().pack(fill="both", expand=True)
            canvas_idx = FigureCanvasTkAgg(fig_idx, master=f_idx); canvas_idx.draw(); canvas_idx.get_tk_widget().pack(fill="both", expand=True)

        ttk.Button(win, text="Cerrar", command=win.destroy).grid(row=2, column=0, sticky="ew", pady=6)

# ---------------------------
# arranque
# ---------------------------
if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()
