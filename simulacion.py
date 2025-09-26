import random
import numpy as np
import matplotlib.pyplot as plt

x0 = [0]*729
rango = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0]

V0 = 1

k = 0
for i in range (0, 81): 
    if (i % 10 == 9): k += 1
    for j in range (0, 9):
        if (i == 0 or (j == 8 and i <= 36)):
            x0[j*81 + i] = V0
        elif (j == 0 or i == 80 or (i >= 58 and j == 8) or (i > 36 and i < 58 and j > 5 and j <= 8) or (i > 10 and i < 20 and j > 0 and j < 3)):
            x0[j*81 + i] = 0
        else:
            x0[j*81 + i] = random.uniform(rango[k], rango[k+1])



""" for j in reversed(range(9)):   # de arriba (NY-1) hacia abajo (0)
    fila = x0[j*81:(j+1)*81]
    fila_str = " ".join(f"{val:4.1f}" for val in fila)
    print(fila_str) """

    # --- visualización ---
NX, NY = 81, 9
matriz = np.array(x0).reshape(NY, NX)

# crear un mapa de calor
plt.figure(figsize=(15, 3))
cmap = plt.cm.viridis  # colormap principal
cmap.set_under('lightgray')  # valores menores al vmin → gris

# plot con imshow
im = plt.imshow(matriz[::-1], cmap=cmap, vmin=0.0001, vmax=1)  
# matriz[::-1] → para que se vea de abajo hacia arriba (como ejes cartesianos)

plt.colorbar(im, label="u valores")
plt.title("Distribución de valores en la malla")
plt.xlabel("i (x)")
plt.ylabel("j (y)")
plt.show()