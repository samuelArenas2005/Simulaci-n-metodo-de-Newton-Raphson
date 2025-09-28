import random
import numpy as np
import matplotlib.pyplot as plt


class Vector:

    def __init__(self,n):
        self.vec = [0]*n
    
    def vecInicial(self,velocidad_init=1):
        V0 = velocidad_init
        k = 0
        rango = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0]
        x0 = self.vec
        for i in range (0, 81): 
            if (i % 10 == 9 and rango[k+1] != 0): k += 1
            for j in range (0, 9):
                if (i == 0 or (j == 8 and i <= 36)):
                    x0[j*81 + i] = V0
                elif (j == 0 or i == 80 or (i >= 58 and j == 8) or (i > 36 and i < 58 and j > 5 and j <= 8) or (i > 10 and i < 20 and j > 0 and j < 3)):
                    x0[j*81 + i] = 0
                elif (i >= 58 and  j > 5 and j <= 8):
                    x0[j*81 + i] = random.uniform(rango[9], rango[10])
                elif (i >= 20 and  j > 0 and j < 3):
                    x0[j*81 + i] = random.uniform(rango[9], rango[10])
                else:
                    x0[j*81 + i] = random.uniform(rango[k], rango[k+1])
                    
                    
                    
                    
                    
    def newVector(self):
        print("Aqui va la logica recursiva")
        
        
                    
    def showInConsole(self,condicion=False):
        x0 = self.vec
        if condicion:
            for j in reversed(range(9)):  
                fila = x0[j*81:(j+1)*81]
                fila_str = " ".join(f"{val:4.1f}" for val in fila)
                print(fila_str)
                
                
    def showPlot(self,condicion=True):
        if condicion:
            x0 = self.vec
                # --- visualización ---
            NX, NY = 81, 9
            matriz = np.array(x0).reshape(NY, NX)

                # crear un mapa de calor
            plt.figure(figsize=(15, 3))
            cmap = plt.cm.viridis  # colormap principal
            cmap.set_under('gray')  # valores menores al vmin → gris
            # valores menores al vmin → gris

            # plot con imshow
            im = plt.imshow(matriz[::-1], cmap=cmap, vmin=0.0001, vmax=1)  
            # matriz[::-1] → para que se vea de abajo hacia arriba (como ejes cartesianos)

            plt.colorbar(im, label="u valores")
            plt.title("Distribución de valores en la malla")
            plt.xlabel("i (x)")
            plt.ylabel("j (y)")
            plt.show()