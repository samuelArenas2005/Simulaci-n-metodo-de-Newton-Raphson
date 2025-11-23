import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline,RegularGridInterpolator

class NewtonRaphson:

    #Valores iniciales de nuestra simulación que permite 
    def __init__(self,Nx, Ny, blocks, V0, v_ij, h,v0ParedSuperior):
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
        self.v0ParedSuperior = v0ParedSuperior
    
    # Método que verifica si un nodo esta dentro de un bloque
    def _is_in_block(self, i, j):
        for block in self.blocks:
            x1, x2, y1, y2 = block
            if x1 <= i <= x2 and y1 <= j <= y2:
                return True
        return False
    
    """ FUNCIONES PRINCIPALES PARA LA SOLUCIÓN DE SISTEMA DE ECUACIONES NO LINEALES MEDIANTE EL MÉTODO DE NEWTON RAPHSON"""
                    
    #Establece los valores del vector inicial
    def vecInicial(self):
        for j in range(self.Ny):
            for i in range(self.Nx):
                index = j * self.Nx + i

                # 1. Comprobar si el nodo está dentro de alguno de los bloques.
                if self._is_in_block(i, j) or i == self.Nx - 1 or j == 0 :
                    self.vec[index] = 0
                
                # 2. Establecer condiciones de frontera si no está en un bloque.
                # Borde izquierdo y superior (entrada de flujo).
                elif i == 0 :
                    self.vec[index] = self.V0
                
                #3. Caso aislado para la pared superior
                elif j == self.Ny - 1:
                    self.vec[index] = self.v0ParedSuperior

                # 4. Si no es frontera ni bloque, es un nodo interior del fluido.
                else:
                    # Se crea un gradiente de velocidad lineal de izquierda a derecha
                    # para tener una suposición inicial más realista.
                    gradiente_lineal = self.V0 * (1 - (i / (self.Nx - 1)))
                    
                    # Se añade una pequeña variación aleatoria para romper la simetría.
                    variacion = self.V0 * random.uniform(-0.05, 0.05)
                    
                    self.vec[index] = gradiente_lineal + variacion

    #Calcula el vector F(uk) de la iteración k 
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

    #Calcula el jacobiano J(uk) de la iteración k 
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
                    
                    
    #Calcula el valor de H mediante el gradienteConjugado como método para la solución de sistema de ecuaciones lineales
    def gradienteConjugado(self, M, e):
        #M es numero de paradas
        j= self.matrixJacobiana
        jt= np.transpose(j)
        
        #Valores iniciales del método
        A= np.matmul(jt,j) #Nuestra nueva matrix A, 
                            #obtenida de multiplicar la transpuerta del jacobiano por el jacobiano para que converga
        b= np.matmul(-jt,self.vectFunction) 
        
        rOld= b
        rNew= 0
        vNew = 0
        vOld= rOld.copy()
        
    
        tk= 0
        xNew= np.zeros(self.Nx*self.Ny)
        sk=0

        #Iteración principal del método del gradiente conjugado
        for i in range(M-1):
                tk= np.dot(rOld,rOld) / np.dot(vOld, np.matmul(A,vOld)) #Calculo de la tasa de aprendizaje
                xNew=xNew + tk*vOld                                     #Calculo del nuevo vector resultante
                rNew=rOld -tk*np.matmul(A,vOld)                         #Calculo del nuevo residuo
                if np.linalg.norm(rNew) < e:                            #Si el residuo es pequeño entonces para las iteraciones
                    print("Convergio por norma el residuo en iteracion:", i+1)
                    break
                sk=np.dot(rNew,rNew)/np.dot(rOld,rOld)                  #Calculo de la nueva tasa de dirección
                vNew=rNew+sk*vOld                                       #Calculo de la nueva direccion de convergencia
                
                rOld=rNew
                vOld=vNew
        return xNew       #Retorno el valor obtenido que sera el H del método de Newthon Raphson 

    #Calculo del nuevo vector en cada iteración de Newthon Raphson utilizando el sistema de ecuaciones lineales para hallar H 
    def newVector(self): 
        h = self.gradienteConjugado(self.Nx*self.Ny,1e-16)
        self.vec = self.vec + h
        
    #Calculo del nuevo vector en cada iteración de Newthon Raphson utilizando la inversa del jacobiano
    def newVectorInversa(self):
        jacobianaInversa = np.linalg.inv(self.matrixJacobiana)
        self.vec = np.subtract(self.vec,np.matmul(jacobianaInversa,self.vectFunction))
    
    
    """COMPROBACIÓN DE LA CONVERGENCIA PARA MÉTODOS ITERATIVOS PARA LA SOLUCIÓN DE SISTEMAS DE ECUACIONES LINEALES"""
    
    #Método que verifica si una matriz A es diagonalmente dominante (jacobin y gauss-seidel)
    def is_diagonally_dominant(self, A):
        for i in range(len(A)):
            sumRow = 0
            a_ii = 0
            for j in range(len(A[i])):
                if(i != j):
                    sumRow += A[i,j]
                else:
                    a_ii = A[i,j]
            if(sumRow > a_ii):
                return False
        return True
    
    #Chequea si la matriz A del sistema de ecuaciones converge, mediante el teorema de convergencia ||(I-Q^(-1)A)|| < 1 (richarson, jacobi, gauss-seidel)
    def check_convergence_theorem(self, metodo):
        Q = np.eye(len(self.matrixJacobiana))
        norm = np.subtract(Q, self.matrixJacobiana)
        #Sacal la nomra 1, la de la suma de colmnas
        return np.linalg.norm(norm, 1) <= 1
    
    #Determina el numero de condicion de nuestra matriz jacobiana
    def cal_numeroCondicion(self):
        jacobiano = self.matrixJacobiana
        jacobianoTranspuesto = jacobiano.T
        jaconew = (jacobianoTranspuesto + jacobiano)/2
        U,S,T = np.linalg.svd(jaconew)

        numeroSingularMax = np.max(S)
        numeroSingularMin = np.min(S)
        numeroCondicion = numeroSingularMax/numeroSingularMin

        return numeroCondicion,numeroSingularMax,numeroSingularMin

    #Verifica si una matriz es simetrica, necesario para la convergencia de los métodos de krilov (Gradiente desce y Gradiente conjugado)
    def is_Simetric(self):
        transpose = np.transpose(self.matrixJacobiana)
        return np.allclose(transpose, self.matrixJacobiana)

    
    """IMPLEMENTACIÓN SPLINES CUBICOS NATURALES"""
    
    def construir_spline_natural(self):
        """
        Construye un spline cúbico natural 1D a partir del vector de velocidades.
        Utiliza splines cúbicos naturales (segunda derivada = 0 en los bordes)
        interpolando primero en dirección x y luego en dirección y.
        
        Retorna:
            spline_object: Objeto con método evaluate(x, y) para evaluar la velocidad
                    en cualquier punto (x, y) del dominio.
        """
        x = np.arange(self.Nx * self.Ny)
        y = self.vec

        cs = CubicSpline(x, y, bc_type='natural')
        
        return cs