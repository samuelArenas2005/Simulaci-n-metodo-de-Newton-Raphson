import random
import numpy as np
import matplotlib.pyplot as plt


class Vector:

    def __init__(self,n):
        self.vec = np.zeros(n)
        self.vectFunction = np.zeros(n)
        self.matrixJacobiana = np.zeros((n,n))
                    

    def vecInicial(self,velocidad_init=1):
        V0 = velocidad_init
        k = 0
        """ rango = [20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 1, 0] """
        rango = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0]
        x0 = self.vec
        for i in range (0, 81): 
            if (i % 10 == 9 and rango[k+1] != 0.1): k += 1
            for j in range (0, 9):
                if (i == 0 or (j == 8 and i < 36)):
                    x0[j*81 + i] = V0
                elif (j == 0 or i == 80 or (i >= 58 and j == 8) or (i >= 36 and i <= 58 and j >= 5 and j <= 8) or (i > 10 and i < 20 and j > 0 and j < 3)):
                    x0[j*81 + i] = 0
                elif (i >= 58 and  j > 5 and j <= 8):
                    x0[j*81 + i] = random.uniform(rango[9], rango[10])
                elif (i >= 20 and  j > 0 and j < 3):
                    x0[j*81 + i] = random.uniform(rango[k+2], rango[k+3])
                else:
                    x0[j*81 + i] = random.uniform(rango[k], rango[k+1])
        

    def vecInicialUniform(self,velocidad_init=1):
        V0 = velocidad_init
        x0 = self.vec
        for i in range (0, 81):
            for j in range (0, 9):
                if (i == 0 or (j == 8 and i < 36)):
                    x0[j*81 + i] = V0
                elif (j == 0 or i == 80 or (i >= 58 and j == 8) or (i >= 36 and i <= 58 and j >= 5 and j <= 8) or (i > 10 and i < 20 and j > 0 and j < 3)):
                    x0[j*81 + i] = 0
                else:
                    x0[j*81 + i] = 1


    def cal_function(self):
        
        val_Vij=0.01
        x0=self.vec

        for i in range (0, 81): 
            for j in range (0, 9):
                if (i == 0 or (j == 8 and i < 36)):
                    self.vectFunction[j*81 + i] = 0
                elif (j == 0 or i == 80 or (i >= 58 and j == 8) or (i >= 36 and i <= 58 and j >= 5 and j <= 8) or (i > 10 and i < 20 and j > 0 and j < 3)):
                    self.vectFunction[j*81 + i] = 0
                else:
                    self.vectFunction[j*81 + i] = 1/4 *( 
                        x0[j*81 + i+1] + x0[j*81 + i-1] +x0[(j+1)*81 + i] +x0[(j-1)*81 + i] 
                        -(5/2* (x0[j*81 + i])* ((x0[j*81 + i+1]) - (x0[j*81 + i-1])))
                        -(5/2*(val_Vij)* ((x0[(j+1)*81 + i]) - (x0[(j-1)*81 + i]))))-x0[j*81 + i]                 


    def showInConsoleFunction(self,condicion=False):
        x0 = self.vectFunction
        if condicion:
            for j in reversed(range(9)):  
                fila = x0[j*81:(j+1)*81]
                fila_str = " ".join(f"{val:4.1f}" for val in fila)
                print(fila_str)                    
                    
                    
                    
    def newVector(self):
        xn = np.subtract(self.vec, np.matmul(self.matrixJacobiana, self.vectFunction))
        
        self.vec = xn

    def showInConsole(self,condicion=False):
        x0 = self.vec
        if condicion:
            for j in reversed(range(9)):  
                fila = x0[j*81:(j+1)*81]
                fila_str = " ".join(f"{val:4.1f}" for val in fila)
                print(fila_str)
                
                
    def showPlot(self,condicion=False):
        if condicion:
            x0 = self.vec
                # --- visualización ---
            NX, NY = 81, 9
            matriz = np.array(x0).reshape(NY, NX)

                # crear un mapa de calor
            plt.figure(figsize=(15, 3))
            cmap = plt.cm.viridis  # colormap principal
            # valores menores al vmin → gris

            # plot con imshow
            im = plt.imshow(matriz[::-1], cmap=cmap, vmin=0.0001, vmax=1)  
            # matriz[::-1] → para que se vea de abajo hacia arriba (como ejes cartesianos)

            plt.colorbar(im, label="u valores")
            plt.title("Distribución de valores en la malla")
            plt.xlabel("i (x)")
            plt.ylabel("j (y)")
            plt.show()

    def cal_jacobiano(self):
        
        v_ij= 0.5

        vec_0=self.vec
        
        self.matrixJacobiana.fill(0.0)
        
        valores_derechos = [ecuacion * 81 - 1  for ecuacion in range(1, 10)]
        valores_bloque_A = [k * 81 + j for j in range(11,20) for k in range(0,3)]
        valores_bloque_B = [k * 81 + j for j in range(36,58) for k in range(5,8)]
        
        for ecuacion in range (729):
            if ((ecuacion % 81 == 0) or (ecuacion >= 0 and ecuacion<=81) or (ecuacion >= 81*8  and ecuacion<=728) or ecuacion in valores_derechos 
                    or ecuacion in valores_bloque_A or ecuacion in valores_bloque_B):
                
                self.matrixJacobiana[ecuacion,ecuacion] = 1
            
            else:
                #aqui van las derivadas parciales, en la posición [i,i] [i,i+1] [i,i-1].... Tal y como lo definimos
                i = ecuacion #Por comodidad visual 
                
                self.matrixJacobiana[ecuacion,i] = (-5/8)*(vec_0[i+1]-vec_0[i-1])-1
                
                self.matrixJacobiana[ecuacion,i+1] = (1/4) - (5/8)*vec_0[i]
                
                self.matrixJacobiana[ecuacion,i-1] = (1/4) + (5/8)*vec_0[i]
                
                self.matrixJacobiana[ecuacion,i+81] = (1/4) - (5/8)*v_ij
                
                self.matrixJacobiana[ecuacion,i-81] = (1/4) + (5/8)*v_ij
                
    def cal_inv_jacobiano(self):
        self.matrixJacobiana =  np.linalg.inv(self.matrixJacobiana)
                    
                    
"""     def showMatrixJacobiana(self, filename="jacobiano.xlsx", show=False):
        if show:
            try:
                mat = self.matrixJacobiana.toarray()
            except:
                mat = self.matrixJacobiana

        df = pd.DataFrame(mat)

        # timestamp para nombre único
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        new_filename = f"jacobiano_{timestamp}.xlsx"

        df.to_excel(new_filename, index=False, header=False)
        print(f"Jacobiano guardado en {new_filename}") """