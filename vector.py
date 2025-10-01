import random
import numpy as np
import matplotlib.pyplot as plt

class Vector:

    def __init__(self,Nx, Ny, blocks, V0, v_ij, h):
        n = Nx * Ny
        self.Nx = Nx
        self.Ny = Ny
        self.blocks = blocks
        self.V0 = V0
        self.v_ij = v_ij
        self.h = h
        self.vec = np.zeros(n)
        self.vectFunction = np.zeros(n)
        self.matrixJacobiana = np.zeros((n,n))
    
    def _is_in_block(self, i, j):
        for block in self.blocks:
            x1, x2, y1, y2 = block
            if x1 <= i <= x2 and y1 <= j <= y2:
                return True
        return False
                    

    def vecInicial(self):
        for j in range(self.Ny):
            for i in range(self.Nx):
                index = j * self.Nx + i

                # 1. Comprobar si el nodo está dentro de alguno de los bloques.
                if self._is_in_block(i, j) or i == self.Nx - 1 or j == 0:
                    self.vec[index] = 0
                
                # 2. Establecer condiciones de frontera si no está en un bloque.
                # Borde izquierdo y superior (entrada de flujo).
                elif i == 0 or j == self.Ny - 1:
                    self.vec[index] = self.V0

                # 3. Si no es frontera ni bloque, es un nodo interior del fluido.
                else:
                    # Se crea un gradiente de velocidad lineal de izquierda a derecha
                    # para tener una suposición inicial más realista.
                    gradiente_lineal = self.V0 * (1 - (i / (self.Nx - 1)))
                    
                    # Se añade una pequeña variación aleatoria para romper la simetría.
                    variacion = self.V0 * random.uniform(-0.05, 0.05)
                    
                    self.vec[index] = gradiente_lineal + variacion

    def cal_function(self):
        # Para mayor claridad, creamos alias para las variables de la clase.
        u = self.vec
        F = self.vectFunction
        Nx = self.Nx
        Ny = self.Ny
        h = self.h # Espaciado de la malla re-escalada

        # Recorremos cada nodo para construir el vector de residuales F(u).
        for j in range(Ny):
            for i in range(Nx):
                index = j * Nx + i

                # 1. Aplicar condiciones segun la posicion del nodo.
                # Para nodos en fronteras con valor fijo o dentro de bloques,
                # el residual es cero, ya que la ecuacion es u_ij = Cte.
                if self._is_in_block(i, j) or i == 0 or j == 0 or i == Nx - 1 or j == Ny - 1:
                    F[index] = 0
                
                # 3. Para nodos interiores, aplicamos la ecuacion discretizada.
                # La ecuacion residual es F(u) = (Ecuacion Discretizada) - u_ij = 0
                else:
                    # Aplicamos la ecuacion residual F(u) = (RHS de Ec. 7) - u_ij en una sola linea.
                    F[index] = (1/4) * (
                        u[index + 1] + u[index - 1] + u[index + Nx] + u[index - Nx]
                        - (h/2) * u[index] * (u[index + 1] - u[index - 1])
                        - (h/2) * self.v_ij * (u[index + Nx] - u[index - Nx])
                    ) - u[index]

    def cal_jacobiano(self):
        u = self.vec
        J = self.matrixJacobiana
        Nx = self.Nx
        Ny = self.Ny
        h_div_8 = self.h / 8.0 # Factor h/8 para h=5

        # Se reinicia la matriz Jacobiana a cero en cada calculo.
        J.fill(0.0)

        # Recorremos cada nodo para construir la fila correspondiente de la matriz.
        for j in range(Ny):
            for i in range(Nx):
                # El indice de la fila de la matriz (la ecuacion k)
                index = j * Nx + i

                # 1. Para nodos en fronteras o bloques, la ecuacion es u_k = Cte.
                # La derivada de u_k con respecto a u_k es 1. El resto son 0.
                if self._is_in_block(i, j) or i == 0 or j == 0 or i == Nx - 1 or j == Ny - 1:
                    J[index, index] = 1.0
                
                # 3. Para nodos interiores, calculamos las 5 derivadas parciales no nulas.
                else:
                    # Las derivadas parciales se calculan a partir de la ecuacion residual:
                    # F_k = 0.25 * (...) - u_k
                    
                    # Derivada respecto a u_k (termino diagonal)
                    J[index, index] = -h_div_8 * (u[index + 1] - u[index - 1]) - 1.0
                    
                    # Derivada respecto a u_{k+1} (vecino derecho)
                    J[index, index + 1] = 0.25 - h_div_8 * u[index]
                    
                    # Derivada respecto a u_{k-1} (vecino izquierdo)
                    J[index, index - 1] = 0.25 + h_div_8 * u[index]
                    
                    # Derivada respecto a u_{k+Nx} (vecino superior)
                    J[index, index + Nx] = 0.25 - h_div_8 * self.v_ij
                    
                    # Derivada respecto a u_{k-Nx} (vecino inferior)
                    J[index, index - Nx] = 0.25 + h_div_8 * self.v_ij

    def cal_inv_jacobiano(self):
        self.matrixJacobiana =  np.linalg.inv(self.matrixJacobiana)

    def newVector(self, alpha=1.0):
        # Calcular la inversa del Jacobiano
        jacobiano_inv = np.linalg.inv(self.matrixJacobiana)
        # x_{n+1} = x_n - alpha * J^(-1) * F(x_n)
        xn = np.subtract(self.vec, alpha * np.matmul(jacobiano_inv, self.vectFunction))
        self.vec = xn
                
                
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